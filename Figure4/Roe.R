############# Roe analysis ################

ziji$plotTissueDist <- function (OR.mtx, k = 2, method.distance = "cosine", do.hclust = T, 
          out.prefix = NULL, OR.max = 3, OR.min = 0, col.rid = "rid", 
          col.ht = circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3)), 
          exp.name = expression(italic(OR)), p.tb = NULL, charSig.tb = NULL, 
          pdf.width = 5.5, pdf.height = 10, ...) 
{
  OR.mtx.tmp.txt <- apply(OR.mtx, 2, function(x) {
    ifelse(x > OR.max, sprintf(">%1.0f", OR.max), ifelse(x < 
                                                           OR.min, sprintf("<%1.0f", OR.min), sprintf("%1.2f", 
                                                                                                      x)))
  })
  rownames(OR.mtx.tmp.txt) <- rownames(OR.mtx)
  if (!is.null(p.tb)) {
    p.mtx <- as.matrix(p.tb[, -c(col.rid), with = F])
    rownames(p.mtx) <- p.tb[[col.rid]]
    p.mtx <- p.mtx[rownames(OR.mtx.tmp.txt), colnames(OR.mtx.tmp.txt), 
                   drop = F]
    p.mtx.txt <- apply(p.mtx, 2, function(x) {
      ifelse(x < 0.05, ifelse(x < 0.01, "**", "*"), "")
    })
    p.mtx.txt[is.na(p.mtx.txt)] <- ""
    show.tmp.txt <- matrix(paste0(OR.mtx.tmp.txt, p.mtx.txt), 
                           nrow = nrow(OR.mtx.tmp.txt))
    rownames(show.tmp.txt) <- rownames(OR.mtx.tmp.txt)
    colnames(show.tmp.txt) <- colnames(OR.mtx.tmp.txt)
  }
  else if (!is.null(charSig.tb)) {
    charSig.mtx <- as.matrix(charSig.tb[, -c(col.rid), with = F])
    rownames(charSig.mtx) <- charSig.tb[[col.rid]]
    charSig.mtx <- charSig.mtx[rownames(OR.mtx), colnames(OR.mtx), 
                               drop = F]
    show.tmp.txt <- charSig.mtx
  }
  else {
    show.tmp.txt <- OR.mtx.tmp.txt
  }
  OR.mtx.tmp <- OR.mtx
  OR.mtx.tmp[OR.mtx.tmp > OR.max] <- OR.max
  OR.mtx.tmp[OR.mtx.tmp < OR.min] <- OR.min
  OR.hclust.row <- run.cutree(OR.mtx.tmp, k = k, method.distance = method.distance, 
                              method.hclust = "ward.D2")
  OR.hclust.row$branch <- dendextend::set(OR.hclust.row$branch, 
                                          "branches_lwd", 2)
  sscVis:::plotMatrix.simple(OR.mtx.tmp, col.ht = col.ht, 
                             out.prefix = sprintf("%s.tissue.dist.rClust.withDend", 
                                                  out.prefix), returnHT = T, show.number = show.tmp.txt,   #  show.number控制热图格子里面 要不要显示数字
                             show.dendrogram = do.hclust, row_dend_width = unit(1.5, "cm"), exp.name = exp.name, 
                             z.hi = OR.max, z.lo = OR.min, mytitle = "Tissue Distribution", 
                             par.heatmap = list(cex.row = 1.5, row_names_gp = grid::gpar(fontsize = 10)), 
                             pdf.width = pdf.width, pdf.height = pdf.height, ...) 
}



library("epitools")
# 计算各个组织的实际观测值
real1 <- table(adata$CellType,adata$Cancer)
expect1 <- expected(table(adata$CellType,adata$Cancer))
roe <- real1/expect1
ziji$plotTissueDist(roe)


# 修改 指定顺序
# adata$Organ <- factor(adata$Organ, levels = c("NOR", "LSIL", "HSIL", "SCC"))
# 修改 指定颜色
real1 <- table(adata$CellType,adata$Cancer)
expect1 <- expected(table(adata$CellType,adata$Cancer))
roe <- real1/expect1
my_colors <- rev(brewer.pal(9, "RdYlBu")) 
col_fun <- colorRamp2(c(0, 1, 3), my_colors[c(1, 5, 9)])
ziji$plotTissueDist( OR.mtx = roe, col.ht = col_fun, exp.name =expression(italic(OR)))


real1 <- table(adata$Cancer,adata$CellType)
expect1 <- expected(table(adata$Cancer,adata$CellType))
roe <- real1/expect1
my_colors <- rev(brewer.pal(9, "RdYlBu")) 
col_fun <- colorRamp2(c(0, 1, 3), my_colors[c(1, 5, 9)])
ziji$plotTissueDist( OR.mtx = roe, col.ht = col_fun, exp.name =expression(italic(OR)))

