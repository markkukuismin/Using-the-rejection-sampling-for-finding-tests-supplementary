
library(tidyverse)
library(gridExtra)

pdist = "norm"

library(tidyverse)
library(gridExtra)

pdist = "norm"
n = c("26", "64", "394")

k = 1

for(i in n){
  
  file_path = paste0("simulations/results/", 
                     pdist, 
                     "_effect_vs_power_", i, "n.txt")
  
  D = read.table(file_path, header = TRUE)
  
  D$Method = ifelse(D$Method == "ar_stat", 
                       "AR test",
                       D$Method)
  
  D$Method = ifelse(D$Method == "t_equalvar", 
                       "t-test",
                       D$Method)
  
  D$Method = as.factor(D$Method)
  
  assign(paste0("Data", k), D)
  
  k = k + 1
  
}

rm(D)

p1 = ggplot(data = Data1, 
           aes(eta, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Difference between means", 
       y = "Statistical power",
       tag = "(A)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"),
        legend.position = "bottom") + 
  guides(shape = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.05, 
             linetype="dashed",
             linewidth = 1)

p2 = ggplot(data = Data2, 
           aes(eta, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Difference between means", 
       y = "Statistical power",
       tag = "(B)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"),
        legend.position = "bottom") + 
  guides(shape = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.05, 
             linetype="dashed",
             linewidth = 1)

p3 = ggplot(data = Data3, 
           aes(eta, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Difference between means", 
       y = "Statistical power",
       tag = "(C)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"),
        legend.position = "bottom") + 
  guides(shape = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.05, 
             linetype="dashed",
             linewidth = 1)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 = p1 + theme(legend.position = "bottom") + ggplot2::guides(shape = guide_legend(nrow = 1))

p1

legend = get_legend(p1)

p1 = p1 + theme(legend.position="none")
p2 = p2 + theme(legend.position="none")
p3 = p3 + theme(legend.position="none")

final_plot = grid.arrange(arrangeGrob(p1, p2, p3,
                                      legend, 
                                      ncol = 3, 
                                      nrow = 2, 
                                      layout_matrix = matrix(c(1, 2, 3,
                                                               4, 4, 4),
                                                             ncol = 3, byrow = T),
                                      widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))
)


fn = paste0("simulations/figures/effect_vs_power_", pdist, ".png")

ggsave(fn,
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 200, 
       units = "mm",
       dpi = 700)
