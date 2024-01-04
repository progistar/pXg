install.packages("rsvg")
library(rsvg)

rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig1.svg", width = 890, file = "/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig1.pdf")
rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig2.svg", width = 1820, file = "/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig2.pdf")
rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig3.svg", width = 1820, file = "/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig3.pdf")
rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig4.svg", width = 1820, file = "/Users/gistar/Documents/MyPapers/pXg/MCP/figures/Fig4.pdf")


rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/2.figures/Fig5.svg", width = 1820,
        file = "/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/2.figures/Fig5.pdf")
dev.off()

rsvg_svg("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/2.figures/Fig5.svg", width = 1820,
         file = "/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/2.figures/Fig5.1.svg")

rsvg_pdf("/Users/gistar/Documents/MyPapers/pXg/MCP/figures/FigS1.svg", width = 1820, file = "/Users/gistar/Documents/MyPapers/pXg/MCP/figures/FigS1.pdf")



setwd("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/2.figures")
rsvg_pdf("revFigs/Figure1.svg", width = 890 ,file = "revFigs/Fig1.pdf")
rsvg_pdf("FigS2.svg", width = 1820 ,file = "FigS2.pdf")
rsvg_pdf("FigS3.svg", width = 1820, file = "FigS3_.pdf")


setwd("/Users/gistar/Documents/MyPapers/pXg/MCP/revision_1/1.RevData/5.Synthetic/UniversalSpectrumViewer")
rsvg_pdf("GEVIGTRW.svg", width = 800, height = 8000,file = "GEVIGTRW.pdf")
