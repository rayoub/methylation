
suppressWarnings(suppressPackageStartupMessages({
	library(shiny)
	library(DT)
	library(here)
	library(bslib)
	library(ggplot2)
	library(forcats)
	library(scales)
}))

source(here::here("R", "files.R"))

# ********************************************************************************
# *** functions
# ********************************************************************************

getBatchIds <- function () {

	dir <- here::here("results", "lab")
	files <- list.files(dir, "*_betas.rds")
	batch_ids <- gsub("_betas.rds","", files)

	return (batch_ids)
}

# ********************************************************************************
# *** ui
# ********************************************************************************

ui <- fluidPage(
	titlePanel("Methylation CNS Tumor Classification"),
	tabsetPanel(
		tabPanel("Import data", 
			br(),
			selectInput("dataset", label = "Dataset", choices = ls("package:datasets")),
			verbatimTextOutput("summary"),
			tableOutput("table")
		),
		tabPanel("View data",
			br(),
			selectInput("selectBatch", label = "Select Batch", choices = getBatchIds()),
			DTOutput("sampleTable"),
			plotOutput("scoresPlot", width = "800px")
		)
	)
)

# ********************************************************************************
# *** server
# ********************************************************************************

server <- function (input, output, session) {

	# *** Import tab

	dataset <- reactive({
		get(input$dataset, "package:datasets")
	})
	
	output$summary <- renderPrint({
		summary(dataset())
	})

	output$table <- renderTable({
		dataset()
	})

	# *** View tab

	r_scores_t <- reactive({
		message("r_scores_t triggered")
		scores <-loadLabData(input$selectBatch, "scores")
		scores_df <- as.data.frame(scores, row.names = rownames(scores))
		scores_t <- tibble::rownames_to_column(scores_df, var = "sample_id")
	})

	r_scores_j <- reactive({
		message("r_scores_j triggered")
		scores_mc <- r_scores_t() |>
			pivot_longer(
				cols = !(sample_id),
				names_to = "mc",
				values_to = "mc_score"
			) |>
			dplyr::group_by(sample_id) |>
			dplyr::slice_max(mc_score) 

		scores_mcf <- r_scores_t() |>
			pivot_longer(
				cols = !(sample_id),
				names_to = "mc",
				values_to = "mc_score"
			) |>
			dplyr::mutate(
				mcf = mcf_lookup(mc),
				.after = mc
			) |>
			dplyr::group_by(sample_id, mcf) |>
			dplyr::summarize(
				mcf_score = sum(mc_score),
				.groups = "drop_last"
			) |>
			dplyr::slice_max(order_by = mcf_score)

		scores_j <- scores_mc |> 
			dplyr::inner_join(scores_mcf, join_by(sample_id)) |>
			dplyr::mutate(
				mc_descr = MC[mc]
			) 
	})

	output$sampleTable <- renderDT({

		datatable(r_scores_j(), 
			colnames = c("Sample ID", "Methylation Class", "MC Score", "Methylation Class Family","MCF Score", "Methylation Class Description"),
			rowname = FALSE,
			selection = list(mode = "single", selected = 1),
			options = list(
					lengthChange = FALSE,
					pageLength = 10
				)
			) |>
			formatPercentage(columns = c("mc_score", "mcf_score"), digits = 2) 

	}, server = FALSE)
	
	output$scoresPlot <- renderPlot({
		s <- input$sampleTable_rows_selected
		if(length(s)) {
			sample_id <- r_scores_j()[s,"sample_id"] |> pull()
			scores_mc <- r_scores_t() |>
				pivot_longer(
					cols = !(sample_id),
					names_to = "mc",
					values_to = "mc_score"
				) |>
				dplyr::group_by(sample_id) |>
				dplyr::filter(
					sample_id == .env$sample_id,
					mc_score >= 0.01
				) |>
				dplyr::slice_max(mc_score, n = 6) |>
				dplyr::mutate(
					mc_descr = MC[mc]
				) 
			ggplot(scores_mc, aes(x = fct_reorder(mc_descr, mc_score, .desc = TRUE), y = mc_score)) +
				geom_col(fill = "steelblue") +
				scale_y_continuous(labels = scales::percent) + 
				scale_x_discrete(labels = label_wrap(width = 10)) + 
				labs(
					x = "Methylation Class",
					y = "Score",
					title = paste("Methylation Class Scores for", sample_id)
				) +
				theme(
					plot.title = element_text(size = 20, face = "bold"), 
					axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 20)),
					axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 20)),
					axis.text = element_text(family = "mono", size = 12, face = "bold"),
					plot.margin = margin(l = 40, b = 40, t = 40)
				)
		}
	})

	#output$scoresPlot <- renderTable({
	#	s <- input$sampleTable_rows_selected
	#	message("selected value:", s)
	#	if(length(s)) {
	#		sample_id <- r_scores_j()[s,"sample_id"] |> pull()
	#		scores_mc <- r_scores_t() |>
	#			pivot_longer(
	#				cols = !(sample_id),
	#				names_to = "mc",
	#				values_to = "mc_score"
	#			) |>
	#			dplyr::group_by(sample_id) |>
	#			dplyr::slice_max(mc_score, n = 6) |>
	#			dplyr::mutate(
	#				mc_descr = MC[mc]
	#			) |>
	#			dplyr::filter(
	#				sample_id == .env$sample_id
	#			)
	#	}
	#}, digits = 3)


#	observeEvent(r_scores(), {
#		message("event triggered")
#		freezeReactiveValue(input, "selectSample")
#		updateSelectInput(inputId = "selectSample", choices = r_scores()$sample_id)
#	})

#	r_scores_max <- reactive({
#		message("r_scores_max trigged")
#		r_scores() |>
#			pivot_longer(
#				cols = !(sample_id),
#				names_to = "mc",
#				values_to = "mc_score"
#			) |>
#			dplyr::group_by(sample_id) |>
#			dplyr::slice_max(mc_score, n = 6) 
#
			#and then filter by input$selectSample
#	})
}

shinyApp(ui, server) 

