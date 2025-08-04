
library(ggplot2)
library(gridExtra)
library(cowplot)
library(GGally)

set.seed(1)

# 2 groups

r = c(-0.99, 0.51, 0.99)

n = 52
p = 2

k = 1:3

G = rep(c(0, 1), each = n)

ID = rep(1:n, 2)

j = 1

pp2 = vector("list", 3)

for(i in k){
  
  Sigma = matrix(r[i], p, p)
  
  diag(Sigma) = 1
  
  X = MASS::mvrnorm(n = n, 
                    mu = rep(0, p),
                    Sigma = Sigma)
  
  
  data = data.frame(X = c(X),
                    G = as.factor(G),
                    ID = as.factor(ID))
  
  d = paste0("data", j)
  
  assign(d, data)
  
  colnames(X) = c("X1", "X2")
  
  X = as.data.frame(X)
  
  pp2[[j]] = GGally::ggpairs(X,
                             upper = list(continuous = GGally::wrap(ggally_cor, stars = F))) +
    labs(tag = paste0("(", LETTERS[j], ")")) +
    theme(plot.tag = element_text(size = 16,
                                  face = "bold"))
  
  j = j + 1
  
}

rm(data)

p1 = ggplot(data1, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[1]], 2)),
       tag = "(A)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")

p2 = ggplot(data2, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[2]], 2)),
       tag = "(B)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")

p3 = ggplot(data3, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[3]], 2)),
       tag = "(C)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")


final_plot = grid.arrange(arrangeGrob(p1, p2, p3,
                                      ncol = 3, 
                                      layout_matrix = matrix(c(1, 2, 3),
                                                             ncol = 3, byrow = T),
                                      widths = c(2.7, 2.7, 2.7), heights = c(2.5))
)

final_plot = plot_grid(
  ggmatrix_gtable(pp2[[1]]),
  ggmatrix_gtable(pp2[[2]]),
  ggmatrix_gtable(pp2[[3]]),
  ncol = 3
)

fn = paste0("simulations/figures/scatter_plotG2.png")

ggsave(fn,
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)

# 3 groups

r = c(-0.49, 0.51, 0.99)

n = 22
p = 3

k = 1:3

G = rep(c(0, 1, 2), each = n)

ID = rep(1:n, 3)

j = 1

pp3 = vector("list", 3)

for(i in k){
  
  Sigma = matrix(r[i], p, p)
  
  diag(Sigma) = 1
  
  X = MASS::mvrnorm(n = n, 
                    mu = rep(0, p),
                    Sigma = Sigma)
  
  data = data.frame(X = c(X),
                    G = as.factor(G),
                    ID = as.factor(ID))
  
  colnames(X) = c("X1", "X2", "X3")
  
  X = as.data.frame(X)
  
  pp3[[j]] = GGally::ggpairs(X, 
                             upper = list(continuous = GGally::wrap(ggally_cor, stars = F))) +
    labs(tag = paste0("(", LETTERS[j], ")")) +
    theme(plot.tag = element_text(size = 16,
                                  face = "bold"))
  
  d = paste0("data", j)
  
  assign(d, data)
  
  j = j + 1
  
}

rm(data)

p1 = ggplot(data1, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[1]], 2)),
       tag = "(A)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")

p2 = ggplot(data2, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[2]], 2)),
       tag = "(B)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")

p3 = ggplot(data3, 
            aes(x = G, y = X, group = ID)) +
  geom_point(aes(shape = G)) +
  geom_line(aes(color = ID)) +
  labs(x = "G", 
       y = "x",
       title = paste0("Correlation between groups = ", round(r[k[3]], 2)),
       tag = "(C)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.tag = element_text(size = 16,
                                face = "bold")) +
  theme(legend.position="none")


final_plot = grid.arrange(arrangeGrob(p1, p2, p3,
                                      ncol = 3, 
                                      layout_matrix = matrix(c(1, 2, 3),
                                                             ncol = 3, byrow = T),
                                      widths = c(2.7, 2.7, 2.7), heights = c(2.5))
)

final_plot = plot_grid(
  ggmatrix_gtable(pp3[[1]]),
  ggmatrix_gtable(pp3[[2]]),
  ggmatrix_gtable(pp3[[3]]),
  ncol = 3
)

fn = paste0("simulations/figures/scatter_plotG3.png")

ggsave(fn,
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)