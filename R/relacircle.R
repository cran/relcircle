relacircle <-
function(table,type1=FALSE,type2=FALSE,line.col=FALSE,pch=1,pch.col="blue",lty=FALSE){
snp<-data.frame(table$rs,table$rs.CHR,table$rs.gene,table$rs.POS)
colnames(snp)<-c("rs","rs.CHR","rs.gene","rs.POS")
snp<-unique(snp)
snp1<-snp[snp$rs.CHR=="X"|snp$rs.CHR=="Y",]
snp2<-snp[snp$rs.CHR!="X"&snp$rs.CHR!="Y",]
snp2[,5]<-as.numeric(snp2[,2])
snp2<-snp2[order(snp2$V5,snp2$rs.POS),]
snp2<-snp2[,-5]
snp1<-snp1[order(snp1$rs.CHR,snp1$rs.POS),]
snp<-rbind(snp2,snp1)

cpg<-data.frame(table$cg,table$cg.CHR,table$cg.gene,table$cg.POS)
colnames(cpg)<-c("cg","cg.CHR","cg.gene","cg.POS")
cpg<-cpg[order(cpg$cg.CHR,cpg$cg.POS),]
cpg<-unique(cpg)
cpg1<-cpg[cpg$cg.CHR=="X"|cpg$cg.CHR=="Y",]
cpg2<-cpg[cpg$cg.CHR!="X"&cpg$cg.CHR!="Y",]
cpg2[,5]<-as.numeric(cpg2[,2])
cpg2<-cpg2[order(cpg2$V5,cpg2$cg.POS),]
cpg2<-cpg2[,-5]
cpg1<-cpg1[order(cpg1$cg.CHR,cpg1$cg.POS),]
cpg<-rbind(cpg2,cpg1)
genome_circle<-function(r,col,lwd,chr_name,chr_length){
  n_chr_name=length(chr_name)
  n_chr_length=length(chr_length)
  if(n_chr_name!=n_chr_length){message('Error,chr# is not match the num of length')
  }else{
    n_chr=length(unique(chr_name))
  }
  gap_angle=0.02*pi
  total_gap_angle=gap_angle*n_chr 
  angle_after_remove_gap=2*pi-total_gap_angle
 if(total_gap_angle+angle_after_remove_gap>2*pi){message("gap is two large!")}
  angle_each_chr= angle_after_remove_gap/n_chr_name 
  chr<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  chr_plot_start<-0
  chr_plot_end<-angle_each_chr
  for(i in 1:n_chr_name){
    n<-which(chr==chr_name[i])
    theta=seq(chr_plot_start,chr_plot_end,by=0.001)
    x=r*cos(theta)
    y=r*sin(theta)
    lines(x,y,col=col[n],lwd=lwd)
    if(i<n_chr_name){
      if(chr_name[i]==chr_name[i+1]){
        chr_plot_start=chr_plot_start+angle_each_chr
        }
      else{
        chr_plot_start=chr_plot_start+angle_each_chr+gap_angle
        }
    }

    if(i<n_chr_name){
      if(chr_name[i]==chr_name[i+1]){
        chr_plot_end=chr_plot_end+angle_each_chr
      }
      else{
        chr_plot_end=chr_plot_end+angle_each_chr+gap_angle
      }
    }
  }
}

inter_line<-function(r_candidate,angle_candidate,r_genome,angle_genome,lwd,col,lty,pch,pch.col){

  angle_diff=angle_genome-angle_candidate
  r_little_circle_candidate=0.1 
  tan_alpha=r_little_circle_candidate/r_candidate 
  alpha=atan(tan_alpha) 
  #
  r_little_circle_genome=0.1
  tan_alpha2=r_little_circle_genome/r_genome
  alpha2=atan(tan_alpha2)

  same_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    if(angle_candidate<=angle_genome){ 
      alpha=alpha
      alpha2=alpha2
      adj_angle=angle_genome-angle_candidate 
      theta=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
      theta2=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
    }else{
      alpha=-alpha
      alpha2=-alpha2
      adj_angle=-1*(angle_genome-angle_candidate)
      theta=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001)
      theta2=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001) 
    }
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }
   left_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    alpha=alpha
    alpha2=alpha2
    adj_angle=angle_genome-angle_candidate 
    theta=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
    theta2=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }

  left_negative_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    alpha=alpha
    alpha2=alpha2
    adj_angle=2*pi-(angle_candidate-angle_genome) 
    theta=seq(pi/2,(adj_angle/2),by=-0.001) 
    theta2=seq(pi/2,(adj_angle/2),by=-0.001) 
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }
 right_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    alpha=-alpha
    alpha2=-alpha2
    adj_angle=-1*(angle_genome-angle_candidate)
    theta=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001) 
    theta2=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001) 
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }
  right_negative_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    alpha=-alpha
    alpha2=-alpha2
    adj_angle=(2*pi-(angle_genome-angle_candidate)) 
    theta=seq(-pi/2,-(adj_angle/2),by=0.001) 
    theta2=seq(-pi/2,-(adj_angle/2),by=0.001)
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }
  opposite_quadrant_line<-function(angle_genome,angle_candidate){
    return_array=c()
    if(angle_genome-angle_candidate>=0){
      if(angle_genome-angle_candidate<=pi){ 
        alpha=alpha
        alpha2=alpha2
        adj_angle=angle_genome-angle_candidate 
        theta=seq(pi/2,(pi/2-adj_angle/2),by=-0.001)
        theta2=seq(pi/2,(pi/2-adj_angle/2),by=-0.001)
      }else{
        alpha=-alpha
        alpha2=-alpha2
        adj_angle=(2*pi-(angle_genome-angle_candidate))
        theta=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001)
        theta2=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001) 
      }
    }else{
      if(angle_candidate-angle_genome<=pi){    
        alpha=-alpha
        alpha2=-alpha2
        adj_angle=angle_candidate-angle_genome 
        theta=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001) 
        theta2=seq(-pi/2,-(pi/2-adj_angle/2),by=0.001)
      }else{
        alpha=alpha
        alpha2=alpha2
        adj_angle=(2*pi-(angle_candidate-angle_genome)) 
        theta=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
        theta2=seq(pi/2,(pi/2-adj_angle/2),by=-0.001) 
      }

    }
    return_list=list(alpha,alpha2,theta,theta2)
    return(return_list)
  }

  if(angle_candidate>=0&&angle_candidate<=pi/2){
    if(angle_genome>=0&&angle_genome<=pi/2){
      return_list=same_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi/2&&angle_genome<=pi){ 
      return_list=left_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi&&angle_genome<=1.5*pi){ 
      return_list=opposite_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }

    if(angle_genome>(1.5*pi)&&angle_genome<2*pi){ 
      return_list=right_negative_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
  }

  if(angle_candidate>pi/2&&angle_candidate<=pi){
    if(angle_genome>=0&&angle_genome<=pi/2){ 
      return_list=right_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi/2&&angle_genome<=pi){ 
      return_list=same_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi&&angle_genome<=1.5*pi){ 
      return_list=left_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }


    if(angle_genome>(1.5*pi)&&angle_genome<2*pi){ 
      return_list=opposite_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }

  }

  if(angle_candidate>pi&&angle_candidate<=1.5*pi){
    if(angle_genome>=0&&angle_genome<=pi/2){ 
      return_list=opposite_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi/2&&angle_genome<=pi){ 
      return_list=right_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi&&angle_genome<=1.5*pi){ 
      return_list=same_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }


    if(angle_genome>(1.5*pi)&&angle_genome<2*pi){ 
      return_list=left_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
  }
  if(angle_candidate>1.5*pi&&angle_candidate<2*pi){
    if(angle_genome>=0&&angle_genome<=pi/2){ 
      return_list=left_negative_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi/2&&angle_genome<=pi){ 
      return_list=opposite_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }
    if(angle_genome>pi&&angle_genome<=1.5*pi){
      return_list=right_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }


    if(angle_genome>(1.5*pi)&&angle_genome<2*pi){ 
      return_list=same_quadrant_line(angle_genome,angle_candidate)
      alpha=return_list[[1]]
      alpha2=return_list[[2]]
      theta=return_list[[3]]
      theta2=return_list[[4]]
    }

  }

   r_alpha=sqrt(r_candidate^2+r_little_circle_candidate^2)
  x_alpha_center=r_alpha*cos(angle_candidate+alpha)
  y_alpha_center=r_alpha*sin(angle_candidate+alpha)
  points(r_candidate*cos(unique(angle_candidate)),r_candidate*sin(unique(angle_candidate)),col=pch.col,pch=pch)
  x_little_cycle1=x_alpha_center+r_little_circle_candidate*cos(angle_candidate-theta)
  y_little_cycle1=y_alpha_center+r_little_circle_candidate*sin(angle_candidate-theta)
  lines(x_little_cycle1,y_little_cycle1,col=col,lwd=lwd)
  r_alpha2=sqrt(r_genome^2+r_little_circle_genome^2)
  x_alpha_center2=r_alpha2*cos(angle_genome-alpha2)
  y_alpha_center2=r_alpha2*sin(angle_genome-alpha2)
  points(r_genome*cos(unique(angle_genome)),r_genome*sin(unique(angle_genome)),col=pch.col,pch=pch)

  x_little_cycle2=x_alpha_center2-r_little_circle_genome*cos(angle_genome-theta2)
  y_little_cycle2=y_alpha_center2-r_little_circle_genome*sin(angle_genome-theta2)
  lines(x_little_cycle2,y_little_cycle2,col=col,lwd=lwd)

  x_little_cycle1_end=x_little_cycle1[length(x_little_cycle1)]
  y_little_cycle1_end=y_little_cycle1[length(y_little_cycle1)]
  r1_end=sqrt(x_little_cycle1_end^2+y_little_cycle1_end^2)
  sin_end_beta1=abs(y_little_cycle1_end/r1_end)
  temp_beta1=asin(sin_end_beta1)
  if(x_little_cycle1_end>=0&&y_little_cycle1_end>=0){beta1=temp_beta1}
  if(x_little_cycle1_end<=0&&y_little_cycle1_end>=0){beta1=pi-temp_beta1}
  if(x_little_cycle1_end<=0&&y_little_cycle1_end<=0){beta1=pi+temp_beta1}
  if(x_little_cycle1_end>=0&&y_little_cycle1_end<=0){beta1=2*pi-temp_beta1}

  x_little_cycle2_end=x_little_cycle2[length(x_little_cycle2)]
  y_little_cycle2_end=y_little_cycle2[length(y_little_cycle2)]
  r2_end=sqrt(x_little_cycle2_end^2+y_little_cycle2_end^2)
  sin_end_beta2=abs(y_little_cycle2_end/r2_end)
  temp_beta2=asin(sin_end_beta2)
  if(x_little_cycle2_end>=0&&y_little_cycle2_end>=0){beta2=temp_beta2}
  if(x_little_cycle2_end<=0&&y_little_cycle2_end>=0){beta2=pi-temp_beta2}
  if(x_little_cycle2_end<=0&&y_little_cycle2_end<=0){beta2=pi+temp_beta2}
  if(x_little_cycle2_end>=0&&y_little_cycle2_end<=0){beta2=2*pi-temp_beta2}

  if(abs(angle_genome-angle_candidate)<=pi){ 
    if(beta1<=beta2){luoxuan_theta_center=seq(beta1,beta2,by=0.001)
    }else{
      luoxuan_theta_center=seq(beta1,beta2,by=-0.001)
    }
  }else{
    if(beta1<=beta2){luoxuan_theta_center=seq(2*pi+beta1,beta2,by=-0.001)
    }else{
      luoxuan_theta_center=seq(beta1,2*pi+beta2,by=0.001)
    }
  }

  r_luoxuan=r1_end
  n_luoxuan_theta_center=length(luoxuan_theta_center)
  delt_luoxuan=(r2_end-r1_end)/n_luoxuan_theta_center
  for(i in 2:n_luoxuan_theta_center){
    r_luoxuan[i]=r_luoxuan[i-1]+delt_luoxuan
  }
  luoxuan_x=r_luoxuan*cos(luoxuan_theta_center)
  luoxuan_y=r_luoxuan*sin(luoxuan_theta_center)
  lines(luoxuan_x,luoxuan_y,col=col,lwd=lwd,lty=lty)

}


genome_circle_point<-function(chr_name,chr,gene){
  n_chr_name=length(chr_name)
  n_chr_length=length(chr)
  if(n_chr_name!=n_chr_length){message('Error,chr# is not match the num of length')
  }else{
    n_chr=length(unique(chr_name))
  }

  gap_angle=0.02*pi
  total_gap_angle=gap_angle*n_chr 
  angle_after_remove_gap=2*pi-total_gap_angle
  if(total_gap_angle+angle_after_remove_gap>2*pi){message("gap is two large!")}
  angle_each_chr= angle_after_remove_gap/n_chr_name 
  chr_plot_start<-c()
  chr_plot_end<-c()
  chr_plot_start[1]=0
  chr_plot_end[1]=angle_each_chr
  point_theta=c()
  for(i in 1:(n_chr_name-1)){

    if(chr_name[i]==chr_name[i+1]){
      chr_plot_start[i+1]=chr_plot_start[i]+angle_each_chr
    }
    else{
      chr_plot_start[i+1]=chr_plot_start[i]+angle_each_chr+gap_angle
    }
    if(i<n_chr_name){
      if(chr_name[i]==chr_name[i+1]){
        chr_plot_end[i+1]=chr_plot_end[i]+angle_each_chr
      }
      else{
        chr_plot_end[i+1]=chr_plot_end[i]+angle_each_chr+gap_angle
      }
    }
    point_theta[i]<-(chr_plot_start[i]+ chr_plot_end[i])/2
  }
  point_theta[n_chr_name]<-(chr_plot_start[n_chr_name]+ chr_plot_end[n_chr_name])/2
  n<-which(chr==gene)
  point_gene<-point_theta[n]
}
plot(x=NULL,y=NULL,xlim=c(-7,7),ylim=c(-7,7),axes = F,xlab ='',ylab='')

chr<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
palette_chr <- palette_chr <-c("#CDE8D1","#D445E5","#C8E89D","#BB6E77","#6FE0A8","#71E664","#CFE751",
                               "#ED8DC7","#BD8CDA","#7341DB","#E5C558","#7071E0","#7B98DE","#74E3E1",
                               "#DDB3A6","#6C9687","#CFD6E6","#7AB5DB","#D948A0","#D172DE","#E58B56",
                               "#E64468","#CDBC85","#DAB5DC")

cpg_point<-c()
for (i in 1:length(table$cg)) {
  cpg_point[i]<-genome_circle_point(cpg$cg.CHR,cpg$cg,table$cg[i])
}
table$cg_theta<-cpg_point
snp_point<-c()
for (i in 1:length(table$rs)) {
  snp_point[i]<-genome_circle_point(snp$rs.CHR,snp$rs,table$rs[i])
}
table$rs_theta<-snp_point
cpg_theta<-data.frame(table$cg,table$cg_theta)
cpg<-merge(cpg,cpg_theta,by.x = colnames(cpg)[1],by.y=colnames(cpg_theta)[1])
colnames(cpg)<-c("cg","cg.CHR","cg.gene","cg.POS","cg.theta")
cpg<-unique(cpg)
cpg1<-cpg[cpg$cg.CHR=="X"|cpg$cg.CHR=="Y",]
cpg2<-cpg[cpg$cg.CHR!="X"&cpg$cg.CHR!="Y",]
cpg2[,6]<-as.numeric(cpg2[,2])
cpg2<-cpg2[order(cpg2[,6],cpg2[,4]),]
cpg2<-cpg2[,-6]
cpg1<-cpg1[order(cpg1$cg.CHR,cpg1$cg.POS),]
cpg<-rbind(cpg2,cpg1)

gene_theta_cpg<-data.frame(table$cg.gene,table$cg_theta)
gene_theta_cpg$table.cg.gene[gene_theta_cpg$table.cg.gene==""] <- NA
gene_theta_cpg<-na.omit(gene_theta_cpg)
gene_theta_cpg<-tapply(gene_theta_cpg$table.cg_theta, INDEX=gene_theta_cpg$table.cg.gene,FUN=mean)
gene_theta_cpg<-as.data.frame(gene_theta_cpg)
gene_theta_cpg[,2]<-row.names(gene_theta_cpg)
colnames(gene_theta_cpg)<-c("theta","gene")

snp_theta<-data.frame(table$rs,table$rs_theta)
snp<-merge(snp,snp_theta,by.x = "rs",by.y="table.rs")
colnames(snp)<-c("rs","rs.CHR","rs.gene","rs.POS","rs.theta")
snp<-unique(snp)
snp1<-snp[snp$rs.CHR=="X"|snp$rs.CHR=="Y",]
snp2<-snp[snp$rs.CHR!="X"&snp$rs.CHR!="Y",]
snp2[,6]<-as.numeric(snp2[,2])
snp2<-snp2[order(snp2[,6],snp2$rs.POS),]
snp2<-snp2[,-6]
snp1<-snp1[order(snp1$rs.CHR,snp1$rs.POS),]
snp<-rbind(snp2,snp1)

gene_theta_snp<-data.frame(table$rs.gene,table$rs_theta)
gene_theta_snp$table.rs.gene[gene_theta_snp$table.rs.gene==""] <- NA
gene_theta_snp<-na.omit(gene_theta_snp)
gene_theta_snp<-tapply(gene_theta_snp$table.rs_theta, INDEX=gene_theta_snp$table.rs.gene,FUN=mean)
gene_theta_snp<-as.data.frame(gene_theta_snp)
gene_theta_snp[,2]<-row.names(gene_theta_snp)
colnames(gene_theta_snp)<-c("theta","gene")


color<-c("#E8A7A9","#CB45E9","#DF7C57","#DBD0DA","#E64A74","#7496E3","#CCE9DB"
         ,"#E48BCD","#CDAEDE","#97C6E7","#D453BD","#D1DC79","#D3BB91","#72E651"
         ,"#C7E7AF","#628494","#634FDD","#70E6BF","#D3E84C","#75E68C","#AD7FE0"
         ,"#E2AF4D","#6FA68C","#9C5F6A","#70E1E6")
r_cpg=3
genome_circle(r_cpg,col = palette_chr,5,cpg$cg.CHR,cpg$cg.POS)

r_snp=1
genome_circle(r_snp,col = palette_chr,3,snp$rs.CHR,snp$rs.POS)

r_snp_line=2
r_cpg_line=3
if(type1[1]!=FALSE&&type2[1]==FALSE){
  table$type<-type1
  type<-unique(type1)
  if(line.col[1]==FALSE){palette<-color[1:length(type)]}else{
    palette<-line.col
  }
  for (i in 1:length(type1)) {
    n<-which(type==type1[i])
    table$col[i]<-palette[n]
  }
  if(lty[1]==FALSE){lty=1}
  for (i in 1:nrow(table)) {
    inter_line(r_snp_line,table$rs_theta[i],r_cpg_line,table$cg_theta[i],1,col = table$col[i],lty=lty,pch=pch,pch.col=pch.col)
  }
  legend("bottomleft",legend = type,lty = lty,col = palette,bty = "n",cex = 0.7,xpd = F,adj = 0.2,text.width = 0.5,ncol = 3,seg.len = 0.5,x.intersp=0.5)
}else
  {
  if(type2[1]!=FALSE&&type1[1]==FALSE){
    table$qtl_type<-type2
    type_qtl<-unique(type2)
    if(line.col[1]==FALSE){palette<-"gray"}else{
      palette<-line.col
    }
    if(lty[1]==FALSE){
      lty<-array(1:length(type_qtl))
      lty<-as.numeric(lty)
    }
    for (i in 1:length(type2)) {
      for (j in 1:length(type_qtl)) {
        if(type2[i]==type_qtl[j]){table$lty[i]<-lty[j]}
      }
    }
    list_type<-list()
    for (i in 1:length(type_qtl)) {
      list_type[[i]]<-table[table$qtl_type==type_qtl[i],]
    }
    for (i in 1:nrow(table)) {
    inter_line(r_snp_line,table$rs_theta[i],r_cpg_line,table$cg_theta[i],1,col = palette,lty=table$lty[i],pch=pch,pch.col=pch.col)
    }
    legend("bottomleft",legend = type_qtl,lty = lty,col =palette,bty = "n",cex = 0.7,xpd = F,adj = 0.2,text.width = 0.5,ncol = 3,seg.len = 0.5,x.intersp=0.5,y.intersp = 0.5)
     }else
      {
        if(type1[1]!=FALSE&&type2[1]!=FALSE){
          table$sig_type<-type1
          type_sig<-unique(type1)
          type_sig<-type_sig[order(type_sig)]
          if(line.col[1]==FALSE){palette<-color[1:length(type_sig)]}else{
            palette<-line.col
          }
          for (i in 1:length(type1)) {
            for (j in 1:length(type_sig)) {
              if(type1[i]==type_sig[j]){table$col[i]<-palette[j]}
            }
          }
          table$qtl_type<-type2
          type_qtl<-unique(type2)
          if(lty[1]==FALSE){
            lty<-array(1:length(type_qtl))
            lty<-as.numeric(lty)
          }
          for (i in 1:length(type2)) {
            for (j in 1:length(type_qtl)) {
              if(type2[i]==type_qtl[j]){table$lty[i]<-lty[j]}
          }
          }

          table$type<-paste(table$sig_type,"-",table$qtl_type,seq="")
          type<-data.frame(table$sig_type,table$qtl_type,table$col,table$lty,table$type)
          type<-type[!duplicated(type$table.type),]
          type<-type[order(type$table.qtl_type,type$table.sig_type),]
          list_type<-list()
          for (i in 1:length(type_qtl)) {
            list_type[[i]]<-type[type$table.qtl_type==type_qtl[i],]
          }
         for (i in 1:nrow(table)) {
            inter_line(r_snp_line,table$rs_theta[i],r_cpg_line,table$cg_theta[i],1,col = table$col[i],lty=table$lty[i],pch=pch,pch.col=pch.col)
          }
          legend("bottomleft",legend = type$table.type,lty = type$table.lty,col =type$table.col,bty = "n",cex = 0.7,xpd = F,adj = 0.2,text.width = 0.5,ncol = 3,seg.len = 0.5,x.intersp=0.5,y.intersp = 0.5)
         
      }else
          {
            if(line.col[1]==FALSE){palette<-"gray"}else{
              palette<-line.col
            }
            if(lty[1]==FALSE){lty=1}
            for (i in 1:nrow(table)) {
              inter_line(r_snp_line,table$rs_theta[i],r_cpg_line,table$cg_theta[i],1,col = palette,lty=lty,pch=pch,pch.col=pch.col)
            }
          }
      }
}


table_rscg<-data.frame(table$cg,table$rs,table$cg_theta,table$rs_theta)
colnames(table_rscg)<-c("cg","rs","cg_theta","rs_theta")
table_cg<-table_rscg[!duplicated(table_rscg$cg),]
table_rs<-table_rscg[!duplicated(table_rscg$rs),]
r_text_cpg=3.5
for (i in 1:nrow(table_cg)) {
  text(x=r_text_cpg*cos(table_cg$cg_theta[i]),y=r_text_cpg*sin(table_cg$cg_theta[i]),table_cg$cg[i],family="sans",cex = 0.7,srt=((table_cg$cg_theta[i]/pi)*180))
}

r_text_snp=1.5
for (i in 1:nrow(table_rs)) {
  text(x=r_text_snp*cos(table_rs$rs_theta[i]),y=r_text_snp*sin(table_rs$rs_theta[i]),table_rs$rs[i],family="sans",cex = 0.6,srt=((table_rs$rs_theta[i]/pi)*180))
}

r_gene_cpg=4.5
for (i in 1:length(gene_theta_cpg[,1])) {
  text(x=r_gene_cpg*cos(gene_theta_cpg$theta[i]),y=r_gene_cpg*sin(gene_theta_cpg$theta[i]),gsub('[|]',"\n",gene_theta_cpg$gene)[i],family="sans",cex = 0.5,srt=(gene_theta_cpg$theta[i]/pi)*180,col="darkblue")
}


r_gene_snp=0.6
for (i in 1:length(gene_theta_snp[,1])) {
  text(x=r_gene_snp*cos(gene_theta_snp$theta[i]),y=r_gene_snp*sin(gene_theta_snp$theta[i]),gsub('[|]',"\n",gene_theta_snp$gene)[i],family="sans",cex = 0.3,srt=(gene_theta_snp$theta[i]/pi)*180,col="darkblue")
}


legend("bottomright",legend = chr,lty = 1,lwd = 15,col = palette_chr,bty = "n",cex = 0.7,xpd = F,adj = 0.2,text.width = 0.1,ncol = 6,x.intersp=0.3,seg.len=0.1,y.intersp=0.5)
}
