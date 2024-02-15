par(las=1,cex=1.25)
defPar <- par(no.readonly = TRUE)

# Setting the colour for the vessel
myPalette = wes_palette("Zissou1", 16, type = "continuous")
vesselCol <- c(myPalette[1],myPalette[16])
vesselBg <- NULL
for(i in 1:2){
  temp <- as.numeric(hex2RGB(vesselCol[i])@coords[1,])
  temp <- rgb(temp[1],temp[2],temp[3],alpha=0.5)
  vesselBg[i] <- temp
}

# Setting pch
myPch <- array(c(0,2,1,15,17,16),dim=c(3,2),dimnames=list(c("vp","vr","T"),c("State 1","State 2")))
