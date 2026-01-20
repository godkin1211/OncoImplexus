library(shiny)
library(shinydashboard)
library(OncoImplexus)
library(BSgenome)
library(DT)
library(ggplot2)

# Increase file upload size limit
options(shiny.maxRequestSize = 200 * 1024^2)

# UI Definition
ui <- dashboardPage(
    dashboardHeader(title = "OncoImplexus Analysis"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Analysis", tabName = "analysis", icon = icon("dna")),
            menuItem("Fusions/Drivers", tabName = "fusions", icon = icon("crosshairs")),
            menuItem("About", tabName = "about", icon = icon("info-circle"))
        ),
        hr(),
        fileInput("sv_file", "Upload SV VCF (.vcf, .vcf.gz)", accept = c(".vcf", ".gz")),
        fileInput("cnv_file", "Upload CNV VCF (.vcf, .vcf.gz)", accept = c(".vcf", ".gz")),
        selectInput("genome", "Genome Build", choices = c("hg38", "hg19")),
        actionButton("run", "Run Analysis", icon = icon("play"), class = "btn-success"),
        hr(),
        helpText("Note: Analysis may take several minutes depending on SV count.")
    ),
    dashboardBody(
        tabItems(
            tabItem(
                tabName = "analysis",
                fluidRow(
                    valueBoxOutput("box_ct", width = 4),
                    valueBoxOutput("box_cp", width = 4),
                    valueBoxOutput("box_cs", width = 4)
                ),
                fluidRow(
                    box(
                        title = "Status", width = 12, status = "info", solidHeader = TRUE,
                        textOutput("status_text")
                    )
                ),
                fluidRow(
                    box(
                        title = "Genome-wide Dashboard", width = 12, status = "primary", solidHeader = TRUE,
                        plotOutput("plot_dashboard", height = "600px")
                    )
                ),
                fluidRow(
                    box(
                        title = "Mechanism Summary", width = 6, status = "warning",
                        plotOutput("plot_mechanisms")
                    ),
                    box(
                        title = "DNA Repair Mechanisms", width = 6, status = "warning",
                        plotOutput("plot_repair")
                    )
                ),
                fluidRow(
                    box(
                        title = "Detailed Event List", width = 12,
                        DTOutput("table_events")
                    )
                )
            ),
            tabItem(
                tabName = "fusions",
                fluidRow(
                    valueBoxOutput("box_drivers", width = 6),
                    valueBoxOutput("box_fusions", width = 6)
                ),
                fluidRow(
                    box(
                        title = "Impacted Driver Genes", width = 12, status = "danger", solidHeader = TRUE,
                        DTOutput("table_drivers")
                    )
                ),
                fluidRow(
                    box(
                        title = "Predicted Gene Fusions", width = 12, status = "primary",
                        DTOutput("table_fusions")
                    )
                )
            ),
            tabItem(
                tabName = "about",
                box(
                    title = "About OncoImplexus", width = 12,
                    h4("Integrated Detection of Complex Chromosomal Rearrangements"),
                    p("OncoImplexus allows for the detection and classification of Chromothripsis, Chromoplexy, and Chromoanasynthesis from WGS data."),
                    p("Version: 1.0.0")
                )
            )
        )
    )
)

# Server Logic
server <- function(input, output, session) {
    # Reactive values to store results
    values <- reactiveValues(
        results = NULL,
        sv_data = NULL,
        cnv_data = NULL
    )

    output$status_text <- renderText("Ready to analyze.")

    observeEvent(input$run, {
        req(input$sv_file, input$cnv_file)

        # Reset
        values$results <- NULL

        withProgress(message = "Analyzing...", value = 0, {
            # 1. Load Data
            output$status_text <- renderText("Step 1/4: Loading SV and CNV data...")
            incProgress(0.1, detail = "Loading VCFs...")

            tryCatch(
                {
                    # Load Genome Library
                    genome_build <- input$genome
                    if (genome_build == "hg19") {
                        library(BSgenome.Hsapiens.UCSC.hg19)
                        genome_obj <- BSgenome.Hsapiens.UCSC.hg19
                    } else {
                        library(BSgenome.Hsapiens.UCSC.hg38)
                        genome_obj <- BSgenome.Hsapiens.UCSC.hg38
                    }

                    # Load Gene Annotation
                    # Assumes gene files are in the same directory or a known location in the container
                    gene_file <- paste0(genome_build, "_genes.rds")
                    if (file.exists(gene_file)) {
                        gene_data <- readRDS(gene_file)
                    } else if (file.exists(file.path("data", gene_file))) {
                        gene_data <- readRDS(file.path("data", gene_file))
                    } else {
                        # Fallback for dev environment or if missing
                        # Ideally this should error out or use a minimal set
                        showNotification(paste("Gene annotation file not found:", gene_file), type = "warning")
                        gene_data <- GRanges()
                    }

                    values$sv_data <- read_sv_vcf(input$sv_file$datapath, genome = input$genome)
                    values$cnv_data <- read_cnv_vcf(input$cnv_file$datapath)

                    # 2. Detection
                    output$status_text <- renderText("Step 2/5: Detecting Chromoanagenesis...")
                    incProgress(0.3, detail = "Running detection algorithms...")

                    res <- detect_chromoanagenesis(
                        SV.sample = values$sv_data,
                        CNV.sample = values$cnv_data,
                        genome = input$genome,
                        gene_granges = gene_data
                    )

                    # 3. Sequence Analysis
                    output$status_text <- renderText("Step 3/5: Analyzing Breakpoint Sequences...")
                    incProgress(0.5, detail = "Analyzing microhomology and repair...")

                    # Need dataframe for this function
                    sv_df <- as(values$sv_data, "data.frame")
                    res$repair <- analyze_breakpoint_sequences(sv_df, genome_obj)

                    # 4. Fusion and Driver Annotation
                    output$status_text <- renderText("Step 4/5: Annotating Genes...")
                    incProgress(0.7, detail = "Predicting fusions and driver impacts...")

                    # Run annotation if function exists
                    if (exists("annotate_chromoanagenesis")) {
                        res$drivers <- annotate_chromoanagenesis(res, gene_data)
                    } else {
                        res$drivers <- list(driver_hits = NULL, all_impacted_genes = NULL)
                    }

                    values$results <- res

                    output$status_text <- renderText("Analysis Complete!")
                    incProgress(1.0, detail = "Done.")
                },
                error = function(e) {
                    showNotification(paste("Error:", e$message), type = "error")
                    output$status_text <- renderText(paste("Error:", e$message))
                }
            )
        })
    })

    # Outputs

    output$box_ct <- renderValueBox({
        if (is.null(values$results) || is.null(values$results$chromothripsis)) {
            val <- "0"
            sub <- "0 High / 0 Low"
        } else {
            n_high <- values$results$chromothripsis$n_high_confidence
            n_low <- values$results$chromothripsis$n_low_confidence
            val <- as.character(n_high + n_low)
            sub <- paste(n_high, "High /", n_low, "Low")
        }
        # Use a div to ensure proper block formatting
        subtitle_tag <- tagList(
            tags$p("Chromothripsis", style = "font-size: 16px; font-weight: bold; margin-bottom: 2px;"),
            tags$p(sub, style = "font-size: 12px; margin-top: 0;")
        )
        valueBox(val, subtitle_tag, icon = icon("bomb"), color = "red")
    })

    output$box_cp <- renderValueBox({
        if (is.null(values$results) || is.null(values$results$chromoplexy)) {
            val <- "0"
            sub <- "0 Likely / 0 Possible"
        } else {
            n_likely <- values$results$chromoplexy$likely_chromoplexy
            n_possible <- values$results$chromoplexy$possible_chromoplexy
            val <- as.character(n_likely + n_possible)
            sub <- paste(n_likely, "Likely /", n_possible, "Possible")
        }
        subtitle_tag <- tagList(
            tags$p("Chromoplexy", style = "font-size: 16px; font-weight: bold; margin-bottom: 2px;"),
            tags$p(sub, style = "font-size: 12px; margin-top: 0;")
        )
        valueBox(val, subtitle_tag, icon = icon("random"), color = "blue")
    })

    output$box_cs <- renderValueBox({
        if (is.null(values$results) || is.null(values$results$chromoanasynthesis)) {
            val <- "0"
            sub <- "0 Likely / 0 Possible"
        } else {
            n_likely <- values$results$chromoanasynthesis$likely_chromoanasynthesis
            n_possible <- values$results$chromoanasynthesis$possible_chromoanasynthesis
            val <- as.character(n_likely + n_possible)
            sub <- paste(n_likely, "Likely /", n_possible, "Possible")
        }
        subtitle_tag <- tagList(
            tags$p("Chromoanasynthesis", style = "font-size: 16px; font-weight: bold; margin-bottom: 2px;"),
            tags$p(sub, style = "font-size: 12px; margin-top: 0;")
        )
        valueBox(val, subtitle_tag, icon = icon("sync"), color = "green")
    })

    output$plot_dashboard <- renderPlot({
        req(values$results)
        plot_genome_dashboard(values$results, values$sv_data, values$cnv_data)
    })

    output$plot_repair <- renderPlot({
        req(values$results$repair)
        plot_repair_mechanisms(values$results$repair)
    })

    output$plot_mechanisms <- renderPlot({
        req(values$results)
        validate(
            need(!is.null(values$results$dominance), "No dominance analysis available"),
            need(!is.null(values$results$dominance$mechanism_proportions), "No mechanism events detected")
        )
        plot_mechanism_dominance(values$results)
    })

    output$table_events <- renderDT({
        req(values$results)

        all_events <- data.frame(
            Mechanism = character(),
            Event_ID = character(),
            Description = character(),
            Classification = character(),
            Confidence = numeric(),
            stringsAsFactors = FALSE
        )

        # 1. Chromothripsis
        if (!is.null(values$results$chromothripsis) && !is.null(values$results$chromothripsis$classification)) {
            ct <- values$results$chromothripsis$classification
            if (nrow(ct) > 0) {
                # Add check for specific columns to avoid errors
                loc_str <- paste0(ct$chrom, ":", ct$start, "-", ct$end)
                conf <- if ("confidence_score" %in% colnames(ct)) ct$confidence_score else NA

                df_ct <- data.frame(
                    Mechanism = "Chromothripsis",
                    Event_ID = paste0("CT_", ct$chrom),
                    Description = loc_str,
                    Classification = ct$classification,
                    Confidence = round(as.numeric(conf), 3),
                    stringsAsFactors = FALSE
                )
                all_events <- rbind(all_events, df_ct)
            }
        }

        # 2. Chromoplexy
        if (!is.null(values$results$chromoplexy) && !is.null(values$results$chromoplexy$summary)) {
            cp <- values$results$chromoplexy$summary
            if (nrow(cp) > 0) {
                # Format ID and Description
                id <- if ("chain_id" %in% colnames(cp)) paste0("CP_", cp$chain_id) else paste0("CP_", 1:nrow(cp))
                desc <- if ("chromosomes_involved" %in% colnames(cp)) as.character(cp$chromosomes_involved) else "Unknown"
                conf <- if ("confidence_score" %in% colnames(cp)) cp$confidence_score else NA

                df_cp <- data.frame(
                    Mechanism = "Chromoplexy",
                    Event_ID = id,
                    Description = paste("Chr:", desc),
                    Classification = cp$classification,
                    Confidence = round(as.numeric(conf), 3),
                    stringsAsFactors = FALSE
                )
                all_events <- rbind(all_events, df_cp)
            }
        }

        # 3. Chromoanasynthesis
        if (!is.null(values$results$chromoanasynthesis) && !is.null(values$results$chromoanasynthesis$summary)) {
            cs <- values$results$chromoanasynthesis$summary
            if (nrow(cs) > 0) {
                # Format ID and Description
                id <- if ("region_id" %in% colnames(cs)) paste0("CS_", cs$region_id) else paste0("CS_", 1:nrow(cs))
                loc_str <- paste0(cs$chrom, ":", cs$start, "-", cs$end)
                conf <- if ("combined_score" %in% colnames(cs)) cs$combined_score else NA

                df_cs <- data.frame(
                    Mechanism = "Chromoanasynthesis",
                    Event_ID = id,
                    Description = loc_str,
                    Classification = cs$classification,
                    Confidence = round(as.numeric(conf), 3),
                    stringsAsFactors = FALSE
                )
                all_events <- rbind(all_events, df_cs)
            }
        }

        if (nrow(all_events) == 0) {
            return(NULL)
        }

        datatable(all_events,
            options = list(scrollX = TRUE, pageLength = 10),
            filter = "top",
            rownames = FALSE
        )
    })

    # Drivers/Fusions Boxes
    output$box_drivers <- renderValueBox({
        val <- if (!is.null(values$results) && !is.null(values$results$drivers$driver_hits)) nrow(values$results$drivers$driver_hits) else 0
        valueBox(val, "Impacted Driver Genes", icon = icon("dna"), color = "red")
    })

    output$box_fusions <- renderValueBox({
        val <- 0
        if (!is.null(values$results) && !is.null(values$results$fusions)) {
            # Count actual fusions (not just disruptions) if preferred, or total
            val <- nrow(values$results$fusions)
        }
        valueBox(val, "Predicted Fusions/Disruptions", icon = icon("crosshairs"), color = "purple")
    })

    # Drivers Table
    output$table_drivers <- renderDT({
        req(values$results)
        if (is.null(values$results$drivers) || is.null(values$results$drivers$driver_hits) || nrow(values$results$drivers$driver_hits) == 0) {
            return(NULL)
        }
        datatable(values$results$drivers$driver_hits,
            options = list(scrollX = TRUE, pageLength = 5),
            fillContainer = FALSE
        )
    })

    # Fusions Table
    output$table_fusions <- renderDT({
        req(values$results)
        if (is.null(values$results$fusions) || nrow(values$results$fusions) == 0) {
            return(NULL)
        }
        # Select key columns
        df <- values$results$fusions
        cols_to_show <- c("fusion_name", "type", "status", "chrom1", "pos1", "chrom2", "pos2", "is_driver_involved")
        show_cols <- intersect(cols_to_show, colnames(df))

        datatable(df[, show_cols, drop = FALSE],
            options = list(scrollX = TRUE, pageLength = 5),
            fillContainer = FALSE
        )
    })
}

shinyApp(ui, server)
