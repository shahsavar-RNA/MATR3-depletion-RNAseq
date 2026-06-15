library(ggplot2)
library(ggforce)

setwd("D:/")

venn_data <- data.frame(
  x = c(-0.3, 0.3),
  y = c(0, 0),
  r = c(0.5, 0.5),
  label = c("U87", "SH-SY5Y")
)

p <- ggplot() +
  geom_circle(data=venn_data,
              aes(x0=x, y0=y, r=r, fill=label),
              alpha=0.6, color="grey40", linewidth=0.8) +
  scale_fill_manual(values=c("U87"="#1a6faf", "SH-SY5Y"="#f5c800")) +
  annotate("text", x=-0.42, y=0, label="784", size=5, fontface="bold", color="white") +
  annotate("text", x=0, y=0, label="79", size=5, fontface="bold", color="white") +
  annotate("text", x=0.42, y=0, label="162", size=5, fontface="bold", color="black") +
  annotate("text", x=-0.45, y=0.62, label="U87", size=4, fontface="bold", color="#1a6faf") +
  annotate("text", x=0.45, y=0.62, label="SH-SY5Y", size=4, fontface="bold", color="#c9a800") +
  coord_fixed() +
  theme_void() +
  theme(legend.position="none",
        plot.background=element_rect(fill="white", color=NA))

ggsave("Venn_splicing.tiff", plot=p, width=10, height=8, units="cm", dpi=300)
