
# Using R version 4.0.3 (2020-10-10)

# https://rich-iannone.github.io/DiagrammeR/graphviz_and_mermaid.html
# https://mikeyharper.uk/flowcharts-in-r-using-diagrammer/
# https://datascienceplus.com/how-to-build-a-simple-flowchart-with-r-diagrammer-package/

#install.packages(c("DiagrammeR", "DiagrammeRsvg","rsvg"))
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(magrittr)

setwd("C:\\Users\\ssteven\\Dropbox\\Deakin\\Chapter_4_Marine_to_terrestrial\\MS_drafts\\figures")

main <- "digraph flowchart {
      # node definitions with substituted label text
      
      graph [rankdir = TB]
      
      node [fontname = Helvetica, shape = rectangle]
      tab1 [label = '@@1'];
      tab2 [label = '@@2'];
      tab3 [label = '@@3'];
      tab4 [label = '@@4'];
      tab5 [label = '@@5'];
      tab6 [label = '@@6'];
      tab7 [label = '@@7']
      tab8 [label = '@@8'];
      tab9 [label = '@@9'];
      tab10 [label = '@@10'];
      tab11 [label = '@@11'];
      tab12 [label = '@@12'];
  
      # edge definitions with the node IDs
      edge [color = red]
      tab1 -> tab2[label = ''];
      tab2 -> tab3[label = ''];
      tab2 -> tab12[label = ''];
      tab3 -> tab4[label = ''];
      tab4 -> tab5[label = ''];
      tab4 -> tab12[label = ''];
      tab5 -> tab6[label = '']; 
      tab6 -> tab7[label = ''];
      tab7 -> tab8[label = ''];
      tab8 -> tab9[label = ''];
      tab9 -> tab10[label = ''];
      tab10 -> tab11[label = ''];
      tab11 -> tab12[label = ''];
      
}

      [1]: '1: Raw outputs'
      [2]: '2: Convert raw outputs to useable files'
      [3]: '3: Create species lists'
      [4]: '4: Connect processed model outputs to species list'
      [5]: '5: Remove burn-in data'
      [6]: '6: Convert model outputs to indicator inputs (subset, aggregate, sample groups)'
      [7]: '7: Calculate indicators'
      [8]: '8: Aggregate indicators and uncertainty measures'
      [9]: '9: Sample time series'
      [10]: '10: Collate final indicator data'
      [11]: '11: Analysis'
      [12]: '12: Visualisation'
    "

plot(main)

grViz(main) %>%
  export_svg %>% charToRaw %>% rsvg_png("main.png")



