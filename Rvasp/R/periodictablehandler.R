#cat(paste0("elements <- do.call(rbind,list(",paste(apply(elements,1,function(x)paste0("c(\"",paste(x,collapse="\",\""),"\")")),collapse=","),"))"))

plot.periodictable <- function(
  highlights=NULL,
  highlightcol="black",
  highlighttexts=NULL,
  underlay=NULL,
  col="grey80",
  col.underlay="red",
  col.Lanthanoid="grey60",
  col.Actinoid="grey90",
  rectmar=0.05,
  xlab="",
  ylab="",
  xaxs="i",
  yaxs="i",
  text.cex=1,
  rows.cex=0.8,
  rows.col="black",
  round=0,
  typ="",
  ...){
  rectmar <- rep(rectmar,length=4)
  
  plot(c(0,18),c(-9.5,0),typ="n",xlab=xlab,ylab=ylab,asp=1,bty="n",xaxs=xaxs,yaxs=yaxs,xaxt="n",yaxt="n",...)
  normal <- periodictable.getremainingelements(highlights)
  rect(2+rectmar[2],-(6-rectmar[1]),3-rectmar[4],-(5+rectmar[[3]]),col= col.Actinoid,border=NA)
  rect(3,-7.5,18,-8.5,col= col.Actinoid,border=NA)
  rect(2+rectmar[2],-(7-rectmar[1]),3-rectmar[4],-(6+rectmar[[3]]),col=col.Lanthanoid,border=NA)  
  rect(3,-8.5,18,-9.5,col=col.Lanthanoid,border=NA)
  if(!is.null(underlay))
  {
    plot.periodictable.addunderlay(underlay,col=col.underlay)
  }
  if(typ!="n")
  {
    
    plot.periodictable.addelements(normal,col=col,rectmar=rectmar,text.cex=text.cex,rows.cex=rows.cex)
    plot.periodictable.addelements(highlights,rows=highlighttexts,col=highlightcol,rectmar=rectmar,text.cex=text.cex,rows.cex=rows.cex,rows.col=rows.col,round=round)
  }
}

plot.periodictable.addelements <- function(elements,
                        text=periodictable.getelementnames(elements),
                        rows=NULL,
                        rectmar=0.05,
                        col="black",
                        fill="white",
                        rows.col="black",
                        rows.cex=0.8,
                        text.cex=1,
                        round=0){
  rectmar <- rep(rectmar,length=4)
  positions <- periodictable.getelementpositions(elements)
  lapply(1:nrow(positions),FUN=function(i)
  {    
    rect(positions[i,1]-1+rectmar[2],-(positions[i,2]-rectmar[1]),positions[i,1]-rectmar[4],-(positions[i,2]-1+rectmar[[3]]),border=col,col=fill)
    offset <- -0.5
    offset2 <- -0.4
    if(!is.null(rows[[i]]))
    {
      offset <- 0.1
      tmpcol <- rep(rows.col,length.out=length(rows[[i]]))
      for(j in 1:length(rows[[i]]))
        {
        offset2 <- offset2 + 0.7*rows.cex
        r <- rows[[i]][j]
        if(length(grep("^[[:digit:].[:space:]]*$", r))>0)
        {
          r <- round(as.numeric(r),round)
        }
        text(positions[i,1]-0.5,-positions[i,2]+0.5,r,pos=1,col=tmpcol[[j]],offset=offset2,cex=rows.cex)
      }
    }
      
    text(positions[i,1]-0.5,-positions[i,2]+0.5,text[i],pos=3,col=col,offset=offset,cex=text.cex)
    })
}

plot.periodictable.addunderlay<-function(underlay,col="red"){
  positions <- periodictable.getelementpositions(underlay)

  lapply(1:nrow(positions),FUN=function(i)
  { 
  # from leftbottom to righttop
   rect(positions[i,1]-1-0.001,-(positions[i,2])-0.001,positions[i,1]+0.001,-(positions[i,2]-1)+0.001,border=NA,col=col,lwd=0)
  })
}

periodictable.getelementnames <- function(element)
{
  elements <- do.call(rbind,list(c("1","H","1s1"),c("2","He","1s2"),c("3","Li","[He] 2s1"),c("4","Be","[He] 2s2"),c("5","B","[He] 2s2 2p1"),c("6","C","[He] 2s2 2p2"),c("7","N","[He] 2s2 2p3"),c("8","O","[He] 2s2 2p4"),c("9","F","[He] 2s2 2p5"),c("10","Ne","[He] 2s2 2p6"),c("11","Na","[Ne] 3s1"),c("12","Mg","[Ne] 3s2"),c("13","Al","[Ne] 3s2 3p1"),c("14","Si","[Ne] 3s2 3p2"),c("15","P","[Ne] 3s2 3p3"),c("16","S","[Ne] 3s2 3p4"),c("17","Cl","[Ne] 3s2 3p5"),c("18","Ar","[Ne] 3s2 3p6"),c("19","K","[Ar] 4s"),c("20","Ca","[Ar] 4s2"),c("21","Sc","[Ar] 3d 4s2"),c("22","Ti","[Ar] 3d2 4s2"),c("23","V","[Ar] 3d3 4s2"),c("24","Cr","[Ar] 3d5 4s"),c("25","Mn","[Ar] 3d5 4s2"),c("26","Fe","[Ar] 3d6 4s2"),c("27","Co","[Ar] 3d7 4s2"),c("28","Ni","[Ar] 3d8 4s2"),c("29","Cu","[Ar] 3d10 4s"),c("30","Zn","[Ar] 3d10 4s2"),c("31","Ga","[Ar] 3d10 4s2 4p"),c("32","Ge","[Ar] 3d10 4s2 4p2"),c("33","As","[Ar] 3d10 4s2 4p3"),c("34","Se","[Ar] 3d10 4s2 4p4"),c("35","Br","[Ar] 3d10 4s2 4p5"),c("36","Kr","[Ar] 3d10 4s2 4p6"),c("37","Rb","[Kr] 5s"),c("38","Sr","[Kr] 5s2"),c("39","Y","[Kr] 4d 5s2"),c("40","Zr","[Kr] 4d2 5s2"),c("41","Nb","[Kr] 4d4 5s"),c("42","Mo","[Kr] 4d5 5s"),c("43","Tc","[Kr] 4d6 5s"),c("44","Ru","[Kr] 4d7 5s"),c("45","Rh","[Kr] 4d8 5s"),c("46","Pd","[Kr] 4d10"),c("47","Ag","[Kr] 4d10 5s"),c("48","Cd","[Kr] 4d10 5s2"),c("49","In","[Kr] 4d10 5s2 5p"),c("50","Sn","[Kr] 4d10 5s2 5p2"),c("51","Sb","[Kr] 4d10 5s2 5p3"),c("52","Te","[Kr] 4d10 5s2 5p4"),c("53","I","[Kr] 4d10 5s2 5p5"),c("54","Xe","[Kr] 4d10 5s2 5p6"),c("55","Cs","[Xe] 6s"),c("56","Ba","[Xe] 6s2"),c("57","La","[Xe] 5d 6s2"),c("58","Ce","[Xe] 4f2 6s2"),c("59","Pr","[Xe] 4f3 6s2"),c("60","Nd","[Xe] 4f4 6s2"),c("61","Pm","[Xe] 4f5 6s2"),c("62","Sm","[Xe] 4f6 6s2"),c("63","Eu","[Xe] 4f7 6s2"),c("64","Gd","[Xe] 4f7 5d 6s2"),c("65","Tb","[Xe] 4f9 6s2"),c("66","Dy","[Xe] 4f10 6s2"),c("67","Ho","[Xe] 4f11 6s2"),c("68","Er","[Xe] 4f12 6s2"),c("69","Tm","[Xe] 4f13 6s2"),c("70","Yb","[Xe] 4f14 6s2"),c("71","Lu","[Xe] 4f14 5d 6s2"),c("72","Hf","[Xe] 4f14 5d2 6s2"),c("73","Ta","[Xe] 4f14 5d3 6s2"),c("74","W","[Xe] 4f14 5d4 6s2"),c("75","Re","[Xe] 4f14 5d5 6s2"),c("76","Os","[Xe] 4f14 5d6 6s2"),c("77","Ir","[Xe] 4f14 5d7 6s2"),c("78","Pt","[Xe] 4f14 5d9 6s"),c("79","Au","[Xe] 4f14 5d10 6s"),c("80","Hg","[Xe] 4f14 5d10 6s2"),c("81","Tl","[Xe] 4f14 5d10 6s2 6p"),c("82","Pb","[Xe] 4f14 5d10 6s2 6p2"),c("83","Bi","[Xe] 4f14 5d10 6s2 6p3"),c("84","Po","[Xe] 4f14 5d10 6s2 6p4"),c("85","At","[Xe] 4f14 5d10 6s2 6p5"),c("86","Rn","[Xe] 4f14 5d10 6s2 6p6"),c("87","Fr","[Rn] 7s"),c("88","Ra","[Rn] 7s2"),c("89","Ac","[Rn] 6d 7s2"),c("90","Th","[Rn] 6d2 7s2"),c("91","Pa","[Rn] 5f2 6d 7s2"),c("92","U","[Rn] 5f3 6d 7s2"),c("93","Np","[Rn] 5f4 6d 7s2"),c("94","Pu","[Rn] 5f6 7s2"),c("95","Am","[Rn] 5f7 7s2"),c("96","Cm","[Rn] 5f7 6d 7s2"),c("97","Bk","[Rn] 5f9 7s2"),c("98","Cf","[Rn] 5f10 7s2"),c("99","Es","[Rn] 5f11 7s2"),c("100","Fm","[Rn] 5f12 7s2"),c("101","Md","[Rn] 5f13 7s2"),c("102","No","[Rn] 5f14 7s2"),c("103","Lr","[Rn] 5f14 6d 7s2"),c("104","Rf","[Rn] 5f14 6d2 7s2"),c("105","Db","[Rn] 5f14 6d3 7s2"),c("106","Sg","[Rn] 5f14 6d4 7s2"),c("107","Bh","[Rn] 5f14 6d5 7s2"),c("108","Hs","[Rn] 5f14 6d6 7s2"),c("109","Mt","[Rn] 5f14 6d7 7s2")))
  selector <- periodictable.getelementselector(element)
  return(elements[selector,2])
}

periodictable.getelementselector <- function(element)
{
  if(is.null(element))
  {
    return(NULL)
  }
  elements <- do.call(rbind,list(c("1","H","1s1"),c("2","He","1s2"),c("3","Li","[He] 2s1"),c("4","Be","[He] 2s2"),c("5","B","[He] 2s2 2p1"),c("6","C","[He] 2s2 2p2"),c("7","N","[He] 2s2 2p3"),c("8","O","[He] 2s2 2p4"),c("9","F","[He] 2s2 2p5"),c("10","Ne","[He] 2s2 2p6"),c("11","Na","[Ne] 3s1"),c("12","Mg","[Ne] 3s2"),c("13","Al","[Ne] 3s2 3p1"),c("14","Si","[Ne] 3s2 3p2"),c("15","P","[Ne] 3s2 3p3"),c("16","S","[Ne] 3s2 3p4"),c("17","Cl","[Ne] 3s2 3p5"),c("18","Ar","[Ne] 3s2 3p6"),c("19","K","[Ar] 4s"),c("20","Ca","[Ar] 4s2"),c("21","Sc","[Ar] 3d 4s2"),c("22","Ti","[Ar] 3d2 4s2"),c("23","V","[Ar] 3d3 4s2"),c("24","Cr","[Ar] 3d5 4s"),c("25","Mn","[Ar] 3d5 4s2"),c("26","Fe","[Ar] 3d6 4s2"),c("27","Co","[Ar] 3d7 4s2"),c("28","Ni","[Ar] 3d8 4s2"),c("29","Cu","[Ar] 3d10 4s"),c("30","Zn","[Ar] 3d10 4s2"),c("31","Ga","[Ar] 3d10 4s2 4p"),c("32","Ge","[Ar] 3d10 4s2 4p2"),c("33","As","[Ar] 3d10 4s2 4p3"),c("34","Se","[Ar] 3d10 4s2 4p4"),c("35","Br","[Ar] 3d10 4s2 4p5"),c("36","Kr","[Ar] 3d10 4s2 4p6"),c("37","Rb","[Kr] 5s"),c("38","Sr","[Kr] 5s2"),c("39","Y","[Kr] 4d 5s2"),c("40","Zr","[Kr] 4d2 5s2"),c("41","Nb","[Kr] 4d4 5s"),c("42","Mo","[Kr] 4d5 5s"),c("43","Tc","[Kr] 4d6 5s"),c("44","Ru","[Kr] 4d7 5s"),c("45","Rh","[Kr] 4d8 5s"),c("46","Pd","[Kr] 4d10"),c("47","Ag","[Kr] 4d10 5s"),c("48","Cd","[Kr] 4d10 5s2"),c("49","In","[Kr] 4d10 5s2 5p"),c("50","Sn","[Kr] 4d10 5s2 5p2"),c("51","Sb","[Kr] 4d10 5s2 5p3"),c("52","Te","[Kr] 4d10 5s2 5p4"),c("53","I","[Kr] 4d10 5s2 5p5"),c("54","Xe","[Kr] 4d10 5s2 5p6"),c("55","Cs","[Xe] 6s"),c("56","Ba","[Xe] 6s2"),c("57","La","[Xe] 5d 6s2"),c("58","Ce","[Xe] 4f2 6s2"),c("59","Pr","[Xe] 4f3 6s2"),c("60","Nd","[Xe] 4f4 6s2"),c("61","Pm","[Xe] 4f5 6s2"),c("62","Sm","[Xe] 4f6 6s2"),c("63","Eu","[Xe] 4f7 6s2"),c("64","Gd","[Xe] 4f7 5d 6s2"),c("65","Tb","[Xe] 4f9 6s2"),c("66","Dy","[Xe] 4f10 6s2"),c("67","Ho","[Xe] 4f11 6s2"),c("68","Er","[Xe] 4f12 6s2"),c("69","Tm","[Xe] 4f13 6s2"),c("70","Yb","[Xe] 4f14 6s2"),c("71","Lu","[Xe] 4f14 5d 6s2"),c("72","Hf","[Xe] 4f14 5d2 6s2"),c("73","Ta","[Xe] 4f14 5d3 6s2"),c("74","W","[Xe] 4f14 5d4 6s2"),c("75","Re","[Xe] 4f14 5d5 6s2"),c("76","Os","[Xe] 4f14 5d6 6s2"),c("77","Ir","[Xe] 4f14 5d7 6s2"),c("78","Pt","[Xe] 4f14 5d9 6s"),c("79","Au","[Xe] 4f14 5d10 6s"),c("80","Hg","[Xe] 4f14 5d10 6s2"),c("81","Tl","[Xe] 4f14 5d10 6s2 6p"),c("82","Pb","[Xe] 4f14 5d10 6s2 6p2"),c("83","Bi","[Xe] 4f14 5d10 6s2 6p3"),c("84","Po","[Xe] 4f14 5d10 6s2 6p4"),c("85","At","[Xe] 4f14 5d10 6s2 6p5"),c("86","Rn","[Xe] 4f14 5d10 6s2 6p6"),c("87","Fr","[Rn] 7s"),c("88","Ra","[Rn] 7s2"),c("89","Ac","[Rn] 6d 7s2"),c("90","Th","[Rn] 6d2 7s2"),c("91","Pa","[Rn] 5f2 6d 7s2"),c("92","U","[Rn] 5f3 6d 7s2"),c("93","Np","[Rn] 5f4 6d 7s2"),c("94","Pu","[Rn] 5f6 7s2"),c("95","Am","[Rn] 5f7 7s2"),c("96","Cm","[Rn] 5f7 6d 7s2"),c("97","Bk","[Rn] 5f9 7s2"),c("98","Cf","[Rn] 5f10 7s2"),c("99","Es","[Rn] 5f11 7s2"),c("100","Fm","[Rn] 5f12 7s2"),c("101","Md","[Rn] 5f13 7s2"),c("102","No","[Rn] 5f14 7s2"),c("103","Lr","[Rn] 5f14 6d 7s2"),c("104","Rf","[Rn] 5f14 6d2 7s2"),c("105","Db","[Rn] 5f14 6d3 7s2"),c("106","Sg","[Rn] 5f14 6d4 7s2"),c("107","Bh","[Rn] 5f14 6d5 7s2"),c("108","Hs","[Rn] 5f14 6d6 7s2"),c("109","Mt","[Rn] 5f14 6d7 7s2")))
  selector <- element
  if(!is.numeric(selector))
    selector <- match(element,elements[,2])
  return(selector)
}

periodictable.getremainingelements<- function(element)
{
  elementscount<-109
  element <- periodictable.getelementselector(element)
  remaining <- 1:elementscount
  if(!is.null(element))
  {
    if(!is.numeric(element))
    {
      element <- which(element)
    }
    remaining <- remaining[-element]
  }
  return(remaining)
}

periodictable.getelementpositions<-function(element)
{
  elementscount<-109
  elements <- do.call(rbind,list(c("1","H","1s1"),c("2","He","1s2"),c("3","Li","[He] 2s1"),c("4","Be","[He] 2s2"),c("5","B","[He] 2s2 2p1"),c("6","C","[He] 2s2 2p2"),c("7","N","[He] 2s2 2p3"),c("8","O","[He] 2s2 2p4"),c("9","F","[He] 2s2 2p5"),c("10","Ne","[He] 2s2 2p6"),c("11","Na","[Ne] 3s1"),c("12","Mg","[Ne] 3s2"),c("13","Al","[Ne] 3s2 3p1"),c("14","Si","[Ne] 3s2 3p2"),c("15","P","[Ne] 3s2 3p3"),c("16","S","[Ne] 3s2 3p4"),c("17","Cl","[Ne] 3s2 3p5"),c("18","Ar","[Ne] 3s2 3p6"),c("19","K","[Ar] 4s"),c("20","Ca","[Ar] 4s2"),c("21","Sc","[Ar] 3d 4s2"),c("22","Ti","[Ar] 3d2 4s2"),c("23","V","[Ar] 3d3 4s2"),c("24","Cr","[Ar] 3d5 4s"),c("25","Mn","[Ar] 3d5 4s2"),c("26","Fe","[Ar] 3d6 4s2"),c("27","Co","[Ar] 3d7 4s2"),c("28","Ni","[Ar] 3d8 4s2"),c("29","Cu","[Ar] 3d10 4s"),c("30","Zn","[Ar] 3d10 4s2"),c("31","Ga","[Ar] 3d10 4s2 4p"),c("32","Ge","[Ar] 3d10 4s2 4p2"),c("33","As","[Ar] 3d10 4s2 4p3"),c("34","Se","[Ar] 3d10 4s2 4p4"),c("35","Br","[Ar] 3d10 4s2 4p5"),c("36","Kr","[Ar] 3d10 4s2 4p6"),c("37","Rb","[Kr] 5s"),c("38","Sr","[Kr] 5s2"),c("39","Y","[Kr] 4d 5s2"),c("40","Zr","[Kr] 4d2 5s2"),c("41","Nb","[Kr] 4d4 5s"),c("42","Mo","[Kr] 4d5 5s"),c("43","Tc","[Kr] 4d6 5s"),c("44","Ru","[Kr] 4d7 5s"),c("45","Rh","[Kr] 4d8 5s"),c("46","Pd","[Kr] 4d10"),c("47","Ag","[Kr] 4d10 5s"),c("48","Cd","[Kr] 4d10 5s2"),c("49","In","[Kr] 4d10 5s2 5p"),c("50","Sn","[Kr] 4d10 5s2 5p2"),c("51","Sb","[Kr] 4d10 5s2 5p3"),c("52","Te","[Kr] 4d10 5s2 5p4"),c("53","I","[Kr] 4d10 5s2 5p5"),c("54","Xe","[Kr] 4d10 5s2 5p6"),c("55","Cs","[Xe] 6s"),c("56","Ba","[Xe] 6s2"),c("57","La","[Xe] 5d 6s2"),c("58","Ce","[Xe] 4f2 6s2"),c("59","Pr","[Xe] 4f3 6s2"),c("60","Nd","[Xe] 4f4 6s2"),c("61","Pm","[Xe] 4f5 6s2"),c("62","Sm","[Xe] 4f6 6s2"),c("63","Eu","[Xe] 4f7 6s2"),c("64","Gd","[Xe] 4f7 5d 6s2"),c("65","Tb","[Xe] 4f9 6s2"),c("66","Dy","[Xe] 4f10 6s2"),c("67","Ho","[Xe] 4f11 6s2"),c("68","Er","[Xe] 4f12 6s2"),c("69","Tm","[Xe] 4f13 6s2"),c("70","Yb","[Xe] 4f14 6s2"),c("71","Lu","[Xe] 4f14 5d 6s2"),c("72","Hf","[Xe] 4f14 5d2 6s2"),c("73","Ta","[Xe] 4f14 5d3 6s2"),c("74","W","[Xe] 4f14 5d4 6s2"),c("75","Re","[Xe] 4f14 5d5 6s2"),c("76","Os","[Xe] 4f14 5d6 6s2"),c("77","Ir","[Xe] 4f14 5d7 6s2"),c("78","Pt","[Xe] 4f14 5d9 6s"),c("79","Au","[Xe] 4f14 5d10 6s"),c("80","Hg","[Xe] 4f14 5d10 6s2"),c("81","Tl","[Xe] 4f14 5d10 6s2 6p"),c("82","Pb","[Xe] 4f14 5d10 6s2 6p2"),c("83","Bi","[Xe] 4f14 5d10 6s2 6p3"),c("84","Po","[Xe] 4f14 5d10 6s2 6p4"),c("85","At","[Xe] 4f14 5d10 6s2 6p5"),c("86","Rn","[Xe] 4f14 5d10 6s2 6p6"),c("87","Fr","[Rn] 7s"),c("88","Ra","[Rn] 7s2"),c("89","Ac","[Rn] 6d 7s2"),c("90","Th","[Rn] 6d2 7s2"),c("91","Pa","[Rn] 5f2 6d 7s2"),c("92","U","[Rn] 5f3 6d 7s2"),c("93","Np","[Rn] 5f4 6d 7s2"),c("94","Pu","[Rn] 5f6 7s2"),c("95","Am","[Rn] 5f7 7s2"),c("96","Cm","[Rn] 5f7 6d 7s2"),c("97","Bk","[Rn] 5f9 7s2"),c("98","Cf","[Rn] 5f10 7s2"),c("99","Es","[Rn] 5f11 7s2"),c("100","Fm","[Rn] 5f12 7s2"),c("101","Md","[Rn] 5f13 7s2"),c("102","No","[Rn] 5f14 7s2"),c("103","Lr","[Rn] 5f14 6d 7s2"),c("104","Rf","[Rn] 5f14 6d2 7s2"),c("105","Db","[Rn] 5f14 6d3 7s2"),c("106","Sg","[Rn] 5f14 6d4 7s2"),c("107","Bh","[Rn] 5f14 6d5 7s2"),c("108","Hs","[Rn] 5f14 6d6 7s2"),c("109","Mt","[Rn] 5f14 6d7 7s2")))
  selector<- periodictable.getelementselector(element)
  data <- lapply(elements[selector,3],function(x){

    s <- strsplit(x," ")[[1]]
    if(length(grep("[",s[1],fixed=T))>0)
    {
      s <- s[-1]
    }
    s <- strsplit(s,"")
    s <- lapply(s,function(x){if(length(x)<3)
                    x[3]<- 1
                              if(length(x)>3)
                              {
                                x[3]<-as.numeric(x[3])*10+as.numeric(x[4])
                                x <- x[1:3]
                              }
                              
                              return(x)})
    s <- do.call(rbind,s)
    row <- as.numeric(s[nrow(s),1])
    if(s[nrow(s),2]=="d")
      row <- row+1
    col <- sum(as.numeric(s[,3]))
    if(row==1&&col>1)
    {
      col <- col+16
    }
    if((row==2||row==3)&&col>2)
    {
      col <- col+10
    }
    if(row==6 || row==7)
    {
      if(col>2 && col<18)
      {
        row <- row+2.5
        col <- col+1
      }
      else{
        if(col>17)
          col <- col-14
      }
    }
    return(c(col,row))
  })
  data <- do.call(rbind,data)
  #data<- cbind(data,elements[selector,2])
  return(data)
}
