library(shiny)
library(iNEXT)
library(DT)

# UI 程式碼
ui <- fluidPage(
  titlePanel("臺灣河川樣站調查"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "選擇檔案"),
      helpText("請確保檔案為CSV檔且包含河川代碼、locality、scientificName以及organismQuantity等columns。"),
      uiOutput("message"),
      uiOutput("alpha_selector"),
      uiOutput("code"),
      uiOutput("T"),
      actionButton("start_analysis", "建議樣站組合分析")
    ),
    
    mainPanel(
      h3("iNEXT分析結果"),
      tabsetPanel(
        type = "tab",
        tabPanel("稀釋外插曲線", plotOutput("p"),verbatimTextOutput("cat")),
        tabPanel("物種覆蓋度之樣站數量建議", tableOutput("summary")),
        tabPanel("建議樣站組合", DT::dataTableOutput("best_sites"), verbatimTextOutput("ca"))
      )
    )
  )
)

# Server 程式碼
server <- function(input, output) {
  observeEvent(input$file, {
    req(input$file)
    
    # 讀取上傳的 CSV 文件
    data <- read.csv(input$file$datapath, header = TRUE)
    
    # 過濾掉河川代碼為 "0" 和一些錯誤代碼的數據
    data = data[data$河川代碼 != "0", ]
    error_river <- c("114", "134", "242", "250")
    data <- data[!(data$河川代碼 %in% error_river), ]
    
    river <- levels(factor(data$河川代碼))
    
    # 定義 river_function
    river_function = function(i){
      data_river = data[data$河川代碼 == river[i], ]
      
      n_sites = length(levels(factor(data_river$locality)))
      n_species = length(levels(factor(data_river$scientificName)))
      
      row_names = levels(factor(data_river$locality))
      col_names = levels(factor(data_river$scientificName))
      
      data.matrix = matrix(NA, nrow = n_sites, ncol = n_species) 
      row.names(data.matrix) = row_names
      colnames(data.matrix) = col_names
      
      for (i in 1:length(row_names)) {
        for (j in 1:length(col_names)) {
          x = data_river[data_river$locality == row_names[i], ]
          y = x[x$scientificName == col_names[j], ]
          z = sum(y$organismQuantity)
          if(z > 0){
            data.matrix[row_names[i], col_names[j]] = 1
          } else {
            data.matrix[row_names[i], col_names[j]] = 0
          }
        }
      }
      
      data.matrix = t(data.matrix)
      
      all_site_num = colSums(data.matrix)
      able_site = which(all_site_num != 0)
      clean_data = data.matrix[, able_site]
      
      input = c(ncol(clean_data), rowSums(clean_data))
      
      result = iNEXT(input, datatype = "incidence_freq")
      return(result)
    }
    
    results = list()
    
    withProgress(message = '檔案上傳中...', value = 0, {
      for (i in 1:length(river)) {
        results[[i]] = river_function(i)
        incProgress(1/length(river), detail = paste("處理河川代碼:", river[i]))
      }  
    })
    
    # 定義 river_analyizes_function 函數
    river_analyizes_function <- function(lo) {
      data_river <- data[data$河川代碼 == river[lo], ]
      
      n_sites <- length(levels(factor(data_river$locality)))
      n_species <- length(levels(factor(data_river$scientificName)))
      
      row_names <- levels(factor(data_river$locality))
      col_names <- levels(factor(data_river$scientificName))
      
      data.matrix <- matrix(NA, nrow = n_sites, ncol = n_species) 
      row.names(data.matrix) <- row_names
      colnames(data.matrix) <- col_names
      
      for (i in 1:length(row_names)) {
        for (j in 1:length(col_names)) {
          x <- data_river[data_river$locality == row_names[i], ]
          y <- x[x$scientificName == col_names[j], ]
          z <- sum(y$organismQuantity)
          if (z > 0) {
            data.matrix[row_names[i], col_names[j]] <- 1
          } else {
            data.matrix[row_names[i], col_names[j]] <- 0
          }
        }
      }
      
      data.matrix <- t(data.matrix)
      
      all_site_num <- colSums(data.matrix)
      able_site <- which(all_site_num != 0)
      clean_data <- data.matrix[, able_site]
      
      return(clean_data)
    }
    
    # 定義 site_collection
    site_collection = function(x, t){
      y = sample(1:ncol(x), t, replace = FALSE)
      site_result = x[, y]
      return(site_result)
    }
    
    # 初始化結果列表
    all_river <- list()
    
    # 遍歷每條河川並進行後續運算
    for (i in 1:length(river)) {
      all_river[[i]] <- river_analyizes_function(i)
    }
    
    # 在這裡進行任何其他後續運算
    output$message <- renderUI({
      "檔案已上傳"
    })
    
    output$alpha_selector <- renderUI({
      if (!is.null(input$file)) {
        numericInput("alpha", "物種覆蓋度:", min = 0, max = 1, step = 0.01, value = 0.95)
      }
    })
    
    # 定義 all_river_site_selection
    all_river_site_selection = function(x, results, alpha){
      river_1 = river_analyizes_function(x)
      
      a = results[[x]][["iNextEst"]][["coverage_based"]]
      
      b = a[a[,4] == "Observed", 2]
      
      c = a[, 2]
      
      if (max(c) >= alpha){
        d = c[c >= alpha]
        e = min(d)
        f = a[a[,2] == e, 3]
        
        if(f >= ncol(river_1)){
          river_code = river[x]
          total_site_num = ncol(river_1)
          advise = "不建議減少樣站"
        } else {
          if(b < alpha){
            river_code = river[x]
            total_site_num = ncol(river_1)
            advise = "不建議減少樣站"
          } else {
            river_code = river[x]
            total_site_num = ncol(river_1)
            advise = paste0("調查", f, "個樣站時，物種覆蓋度可達", round(e * 100,2), "%")
          }
        }
      } else {
        river_code = river[x]
        total_site_num = ncol(river_1)
        advise = "不建議減少樣站"
      }
      
      return(c(river_code,total_site_num,advise))
      
    }
    
    output$summary <- renderTable({
      req(input$alpha)
      x = matrix(NA,nrow = length(river),ncol = 3)
      colnames(x) = c("河川代碼","總樣站數","建議")
      x = as.data.frame(x)
      for (i in 1:length(results)) {
        y = all_river_site_selection(i,results,input$alpha)
        x[i,1:3] = y[1:3]
      }
      x
    })
    
    output$code <- renderUI({
      if (!is.null(input$alpha)) {
        selectInput("c", "請選擇河川代碼:", river)
      }
    })
    
    output$p = renderPlot({
      if (!is.null(input$c)) {
        a <- as.numeric(input$c)
        b = which(river == a)
        if (!is.na(b)) {
          plot(results[[b]])
        }
      } 
    })
    
    output$T = renderUI({
      if (!is.na(input$c)){
        a <- as.numeric(input$c)
        b = which(river == a)
        river_subset = subset(results[[b]]$iNextEst$coverage_based, Method == "Rarefaction")
        if(!is.na(b)){
          selectInput("t", "請選擇樣站數量:", river_subset[, 3])
        }
      }
    })
    
    output$cat <- renderText({
      river_num <- which(river == input$c)
      co <- results[[river_num]]$iNextEst$coverage_based$t
      paste(
        "物種覆蓋度:", round(results[[river_num]]$iNextEst$coverage_based[co == input$t, 2] * 100, 2), "%")
    })
    
    observeEvent(input$start_analysis, {
      req(input$c, input$t)
      
      # 定義進度條
      withProgress(message = '組合樣站中...', value = 0, {
        river_num <- which(river == input$c)
        
        site_selection_result <- function(x, t, n = 1000000){
          total_site_selection <- list()
          
          for (i in 1:n) {
            total_site_selection[[i]] <- site_collection(x, t)
            if (i %% 10000 == 0) {  
              incProgress(10000 / n, detail = paste("進度:", round(i / n * 100, 2), "%"))
            }
          }
          
          return(total_site_selection)
        }
        
        a <- site_selection_result(all_river[[river_num]], input$t)
        
        total <- sum(all_river[[river_num]])
        
        site_1 <- c()
        
        for (i in 1:nrow(all_river[[river_num]])) {
          site_1[i] <- sum(all_river[[river_num]][i, ])
        }
        
        true_ratio <- site_1 / total
        
        bias_function <- function(x, true_ratio){
          total <- sum(x)
          
          site_1 <- c()
          
          for (i in 1:nrow(x)) {
            site_1[i] <- sum(x[i, ])
          }
          
          ratio <- site_1 / total
          
          bias <- mean(abs(true_ratio - ratio))
          
          return(bias)
        }
        
        k <- c()
        
        for (i in 1:length(a)) {
          k[i] <- bias_function(a[[i]], true_ratio)
        }
        
        w <- which(results[[river_num]]$iNextEst$coverage_based[, 3] == input$t)
        
        n <- as.numeric(results[[river_num]]$iNextEst$coverage_based[w, 6])
        
        co <- results[[river_num]]$iNextEst$coverage_based$t
        
        best_sites <- colnames(a[[which.min(k)]])
        
        output$best_sites <- DT::renderDataTable({
          data.frame("建議樣站" = best_sites)
        }, options = list(pageLength = 10))
      })
      
      output$ca <- renderText({
        paste("以上為建議樣站組合，物種組成差異:", round(min(k) * 100, 2), "%，",
              "物種數:", round(n))
      })
    })
  })
}

# 啟動應用
shinyApp(ui = ui, server = server)