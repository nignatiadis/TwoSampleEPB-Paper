library(RColorBrewer)

Full_FDR_Data = read.csv('./Full Data/Diffuse_FDR_Data.csv')
Full_Power_Data = read.csv('./Full Data/Diffuse_Power_Data.csv')

method_colors = c(
  EV_test = "#D55E00",   # orange (plus)
  Welch  = "#E69F00",    # gold (triangle)
  B_F    = "#CC79A7",    # magenta (diamond)
  VREPB  = "#0072B2",    # blue (square)
  DVEPB  = "#009E73"     # green (circle)
)

method_pch = c(
  EV_test = 3,     # plus sign
  Welch   = 17,    # triangle
  B_F     = 18,    # diamond
  VREPB   = 15,    # square
  DVEPB   = 16     # circle
)

FDR_plotter = function(a, Full_FDR_Data) {
  
  dat = Full_FDR_Data[Full_FDR_Data$a == a, ]
  
  plot(dat$n2, dat$VREPB, type = 'o', pch = method_pch["VREPB"], col = method_colors["VREPB"], 
       lwd = 3, cex = 1.5, xlab = expression(K[B]), ylab = "FDR",
       ylim = c(0, max(dat[, c("VREPB", "DVEPB", "Welch", "B_F", "EV_test")], 0.13)),
       cex.lab = 1.5, cex.axis = 1.5, bty = "n")
  lines(dat$n2, dat$DVEPB, type = 'o', pch = method_pch["DVEPB"], col = method_colors["DVEPB"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$Welch, type = 'o', pch = method_pch["Welch"], col = method_colors["Welch"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$B_F, type = 'o', pch = method_pch["B_F"], col = method_colors["B_F"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$EV_test, type = 'o', pch = method_pch["EV_test"], col = method_colors["EV_test"], lwd = 3, cex = 1.5)
  abline(h = 0.1, lwd = 2, col = "gray40", lty = 2)
}

Power_plotter = function(a, Full_Power_Data) {
  
  dat = Full_Power_Data[Full_Power_Data$a == a, ]
  
  plot(dat$n2, dat$VREPB, type = 'o', pch = method_pch["VREPB"], col = method_colors["VREPB"],
       lwd = 3, cex = 1.5, xlab = expression(K[B]), ylab = "Power",
       ylim = c(0, max(dat[, c("VREPB", "DVEPB", "Welch", "B_F", "EV_test")])), 
       cex.lab = 1.5, cex.axis = 1.5, bty = "n")
  lines(dat$n2, dat$DVEPB, type = 'o', pch = method_pch["DVEPB"], col = method_colors["DVEPB"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$Welch, type = 'o', pch = method_pch["Welch"], col = method_colors["Welch"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$B_F, type = 'o', pch = method_pch["B_F"], col = method_colors["B_F"], lwd = 3, cex = 1.5)
  lines(dat$n2, dat$EV_test, type = 'o', pch = method_pch["EV_test"], col = method_colors["EV_test"], lwd = 3, cex = 1.5)
}

png("./Graph/Pub Combined Diffuse labeled.png", width = 3200, height = 2400, res = 300, family = "sans")
par(mfrow = c(3, 2), mar = c(4,4,2,2), oma = c(4, 11, 2, 13), mgp = c(2.8, 1, 0))

# First row (a = 5)
FDR_plotter(5, Full_FDR_Data)
Power_plotter(5, Full_Power_Data)

# First row (a = 10)
FDR_plotter(10, Full_FDR_Data)
Power_plotter(10, Full_Power_Data)

# First row (a = 15)
FDR_plotter(15, Full_FDR_Data)
Power_plotter(15, Full_Power_Data)

# Add row labels on the left, centered vertically for each row
par(xpd = NA)
mtext(expression(bold("a) ") * a == 5),
      side = 2, outer = TRUE, line = 2, at = 0.95, cex = 1.5, las = 1)
mtext(expression(bold("b) ") * a == 10),
      side = 2, outer = TRUE, line = 2, at = 0.62, cex = 1.5, las = 1)
mtext(expression(bold("c) ") * a == 15),
      side = 2, outer = TRUE, line = 2, at = 0.29, cex = 1.5, las = 1)

# Add legend
legend("topright", inset = c(-0.57, -3.15),
       legend = c("EV-test", "Welch", "B-F", "VREPB", "DVEPB"),
       col = method_colors, pch = method_pch, lwd = 3, pt.cex = 2,
       cex = 1.5, y.intersp = 1.3, bty = "n", xpd = NA)

dev.off()