
library(tidyverse)
library(gridExtra)

pdist = "norm" # norm, laplace, t
d = c("04", "06", "08")

k = 1

for(i in d){
  
  file_path = paste0("simulations/results/", 
                     pdist, 
                     "_52n_",
                     i,
                     "d.txt")
  
  D = read.table(file_path, 
                    header = TRUE)
  
  D$Method = ifelse(D$Method == "ar_stat", 
                    "AR stat.",
                    D$Method)
  
  D$Method = ifelse(D$Method == "t_paired", 
                    "Paired t-test",
                    D$Method)
  
  D$Method = ifelse(D$Method == "t_equalvar", 
                    "t-test",
                    D$Method)
  
  D$Method = as.factor(D$Method)

  if(pdist == "norm"){
    
    n = 52
    
    r = seq(-0.99, 0.99, length.out = 50)
    
    r = c(r, 0)
    
    r = sort(r)
    
    b = 2*(1 - r)/n
    
    cdot = sqrt(b)*qnorm(0.95)
    
    eta = as.numeric(i)/10
    
    ef = (1/sqrt(b))*(cdot - eta)
    
    lr_power = 1 - pnorm(ef)
    
    Data_temp = data.frame(Method = "LR",
                           power = lr_power,
                           error = 0.05,
                           r = r)
    
    Data_temp$Method = as.factor(Data_temp$Method)
    
    D = rbind(D, Data_temp)
    
  }
  
  assign(paste0("Data", k), D)
    
  k = k + 1
  
}

rm(D)

p1 = ggplot(data = Data1, 
            aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power",
       title = paste0("Difference between pop. means = ", as.numeric(d[1])/10),
       tag = "(A)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p2 = ggplot(data = Data2, 
            aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power",
       title = paste0("Difference between pop. means = ", as.numeric(d[2])/10),
       tag = "(B)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p3 = ggplot(data = Data3, 
            aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power",
       title = paste0("Difference between pop. means = ", as.numeric(d[3])/10),
       tag = "(C)") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p4 = ggplot(data = Data1, 
            aes(r, error, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Type I error",
       title = "Significance level = 0.05",
       tag = "(D)") +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.05,
             linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

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
p4 = p4 + theme(legend.position="none")

final_plot = grid.arrange(arrangeGrob(p1, p2,
                                      p3, p4,
                                      legend, 
                                      ncol = 2, 
                                      nrow = 3, 
                                      layout_matrix = matrix(c(1, 2,
                                                               3, 4,
                                                               5, 5),
                                                             ncol = 2, byrow = T),
                                      widths = c(2.7, 2.7), heights = c(2.5, 2.5, 0.2))
)

fn = paste0("simulations/figures/power_", pdist, ".png")

ggsave(fn,
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)
