library(ggplot2)

setwd("D:/")

data <- data.frame(
  Event = rep(c("SE","A3SS","A5SS","MXE","RI"), 2),
  Count = c(208, 13, 13, 17, 5, 679, 79, 65, 79, 59),
  CellLine = c(rep("SH-SY5Y",5), rep("U87",5))
)

data$Event <- factor(data$Event, levels=c("SE","A3SS","A5SS","MXE","RI"))
data$CellLine <- factor(data$CellLine, levels=c("SH-SY5Y","U87"))

p <- ggplot(data, aes(x=Event, y=Count, fill=CellLine)) +
  geom_bar(stat="identity", position=position_dodge(width=0.65),
           width=0.6, color="black", linewidth=0.3) +
  scale_fill_manual(values=c("SH-SY5Y"="#1a6faf", "U87"="#f5c800"),
                    labels=c("SH-SY5Y (neuroblastoma)","U87 (glioblastoma)")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,750), breaks=seq(0,700,100)) +
  labs(x="Alternative splicing event type",
       y="Significant events (FDR < 0.05)", fill=NULL) +
  theme_classic(base_size=11) +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    panel.background=element_rect(fill="white", color=NA),
    axis.line=element_line(color="black", linewidth=0.5),
    axis.title.x=element_text(size=9, color="black", margin=margin(t=8)),
    axis.title.y=element_text(size=9, color="black", margin=margin(r=8)),
    axis.text=element_text(size=9, color="black"),
    legend.position="top",
    legend.key.size=unit(0.4,"cm"),
    legend.text=element_text(size=8),
    panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
    plot.margin=margin(15,20,10,15)
  )

ggsave("Splicing_barplot.tiff", plot=p, width=13, height=9, units="cm", dpi=300)
