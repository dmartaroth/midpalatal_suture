# ## ######################################## ## #
#                     THEMES                     #
# ## ######################################## ## #
my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A", "thistle1","#CB4335","#C7CC8F")


pastel_palette <- c("#7EBDC2", "#D99B82", "#C7CC8F", "#B7A4DB", "#FFB5B5", "#8FC1A9", "#FFD966", "#B2A39E", "#C9ADA7", "#A5B9C4",
                    "#E6B0C2", "#7F8C8D", "#FADBD8", "#ABEBC6", "#D5DBDB", "#F5CBA7", "#E59866", "#641E16", "#F8C471", "#D35400",
                    "#2E4053", "#6C3483", "#2980B9", "#7D3C98", "#F4D03F", "#1F618D", "#6E2C00", "#B3B6B7", "#154360", "#FAD7A0",
                    "#9A7D0A", "#873600", "#DC7633", "#4A235A", "#424949", "#8E44AD", "#1B4F72", "#CB4335", "#76448A", "#2E86C1",
                    "#F1C40F", "#F1948A", "thistle3")

visiumcolors <- c("0" = "#4F94CD","1" = "red","2" = "orchid4","3" = "#FFD700",
                  "4" = "navy","5" = "orange","6" = "lavender","7" = "orchid",
                  "8" = "#AB82FF","9" = "#008B00","10" = "#A52A2A","11" = "#00FFFF")


# Custom dotplot theme
custom_dotplot_theme <- function() {
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 8), 
        legend.text = element_text(size = 8),  
        legend.title = element_text(size = 8)  
  )
}

