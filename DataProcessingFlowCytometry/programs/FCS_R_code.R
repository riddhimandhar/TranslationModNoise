args=commandArgs();

library(flowCore);
library(flowViz);
library(MASS);
library(aspace);
library(graphics);
library(colorspace);

filename=args[8];
outfilename=args[9];
dirname=args[10];
imgoutname=args[11];
txtoutname=args[12];

FRAC=0.5;
N=20;
YLIM1=150;
YLIM2=150;

sink(outfilename,append=TRUE,split=FALSE);

frames <- lapply(dir(path=dirname,filename,full.names=TRUE, all.file=FALSE,recursive=FALSE,ignore.case=FALSE,include.dirs=FALSE,no..=FALSE),read.FCS)
fs<-as(frames,"flowSet");
##plot(transform(fs[[1]], `FSC-A`=log(`FSC-A`), `SSC-A`=log(`SSC-A`)),c("FSC-A","SSC-A"))
#3plot(exprs(fs$V1[,"FSC-A"]),exprs(fs$V1[,"SSC-A"]),log="xy")
##dev.off()
qu1=quantile(exprs(fs$V1[,"FSC-A"]),probs=seq(0,1,0.1));
qu2=quantile(exprs(fs$V1[,"SSC-A"]),probs=seq(0,1,0.1));
YLIM1=qu1["90%"]; 
YLIM2=qu2["90%"]; 

##N=max(YLIM1,YLIM2)/150*20;

N=100;

dist.interval1=(YLIM1/(N-1));
dist.interval2=(YLIM2/(N-1));

f <- kde2d((exprs(fs$V1[,"FSC-A"])), (exprs(fs$V1[,"SSC-A"])), n = N,lims=c(0,YLIM1,0,YLIM2));

cov.matrix <- function (a, b, angle) {
theta <- angle * (pi/180);

c1 <- ((cos(theta)^2)/a^2) + ((sin(theta)^2)/b^2);
c2 <- sin(theta) * cos(theta) * ((1/a^2) - (1/b^2));
c3 <- ((sin(theta)^2)/a^2) + ((cos(theta)^2)/b^2);

m1 <- matrix(c(c1, c2, c2, c3), byrow=TRUE, ncol=2);
m2 <- solve(m1);
m2;
}

sorted.fz=sort(f$z,decreasing=TRUE);
cnt=N*N;
totsum=sum(sorted.fz[1:cnt]);

fraccount=0.05;

while(fraccount<=FRAC)
{

print(fraccount);

i=1;
runsum=0;
lim=0;
while(i <= cnt)
{
  if(lim==0) { runsum=runsum+sorted.fz[i]; }  
  if(lim==0 && runsum>(fraccount*totsum)) { lim=i; }
  i=i+1;
}

##if(lim==2) { lim=3; }

xpos=0; ypos=0;
i=1;
while(i<=lim)
{
   maxp=which(f$z==sorted.fz[i],arr.ind=TRUE);
   xpos=xpos+f$x[maxp[1,1]]+(dist.interval1/2)
   ypos=ypos+f$y[maxp[1,2]]+(dist.interval2/2)
   i=i+1
}

xpos=xpos/lim
ypos=ypos/lim

ptxlist=array()
ptylist=array()

i=1
num=1
while(i<=lim)
{
  p=which(f$z==sorted.fz[i],arr.ind=TRUE)
  ptxlist[num]=f$x[p[1,1]]
  ptylist[num]=f$y[p[1,2]]
  ptxlist[num+1]=f$x[p[1,1]]+dist.interval1
  ptylist[num+1]=f$y[p[1,2]]
  ptxlist[num+2]=f$x[p[1,1]]
  ptylist[num+2]=f$y[p[1,2]]+dist.interval2
  ptxlist[num+3]=f$x[p[1,1]]+dist.interval1
  ptylist[num+3]=f$y[p[1,2]]+dist.interval2
  num=num+4
  i=i+1
}

i=1
maxdist=0
maxdx=0; maxdy=0;
while(i<num)
{
   x=ptxlist[i]
   y=ptylist[i]
   dist=sqrt((xpos-x)^2+(ypos-y)^2)
   if(dist>maxdist) { maxdist=dist; maxdx=x; maxdy=y; }
   i=i+1
}

angle=atan_d((maxdy-ypos)/(maxdx-xpos))

i=1
mindist=0; min.angle=0;
mindx=0; mindy=0;
while(i<num)
{
   x=ptxlist[i]
   y=ptylist[i]
   dist=sqrt((xpos-x)^2+(ypos-y)^2)
   calc.angle=atan_d((y-ypos)/(x-xpos))
   if(dist>mindist && ((calc.angle>=angle+90-45 && calc.angle<=angle+90+45) || (calc.angle>= angle-90-45 && calc.angle<=angle-90+45))) { mindist=dist; min.angle=calc.angle; mindx=x; mindy=y; }
   i=i+1
}

if(abs(angle)>90) { angle=180-(abs(angle)) }
if(abs(calc.angle)>90 ) { calc.angle=180-(abs(calc.angle)) }

lengthfrac=0.75;
while(lengthfrac<=1)
{

majaxis=maxdist*lengthfrac;

minaxis=mindist*cos((90-abs(min.angle)-abs(angle))*pi/180)*lengthfrac;
if(minaxis==0) { minaxis=min(dist.interval1,dist.interval2)*lengthfrac; }


##max(f$z)
cov=cov.matrix(majaxis,minaxis,angle)
colnames(cov) <- c("FSC-A", "SSC-A")
rownames(cov) <- c("FSC-A", "SSC-A")

mean <- c("FSC-A"=xpos, "SSC-A"=ypos)
eg <- ellipsoidGate(filterId= "myEllipsoidGate", .gate=cov, mean=mean)
fres <- filter(fs, eg);
fres;
summary(fres);

filtdat=Subset(fs, fres);

imgfilename=paste(imgoutname,fraccount,lengthfrac,sep="_");
#imgfilename=paste(imgfilename,"jpg",sep=".");
imgfilename=paste(imgfilename,"pdf",sep=".");

#jpeg(imgfilename,width=800,height=600,quality=100);
pdf(imgfilename);
#contour(f, xlab = "FSC-A", ylab = "SSC-A",cex.lab=1.3,cex.axis=1.3,font.lab=2, font.axis=2,family="serif");
#points(exprs(filtdat$V1[,"FSC-A"]),exprs(filtdat$V1[,"SSC-A"]),col="red");
zlim=range(f$z,finite=TRUE);
lev=length(pretty(zlim,N))-1
color.palette.terrain=rev(terrain.colors(lev));
#color.palette.seq=rev(sequential_hcl(n=lev, h = 260, c = c(80, 0), l = c(30, 90), power = 1.5));
#color.palette.heat=rev(heat_hcl(n=lev, h = c(0, 90), c = c(100, 30), l = c(50, 90), power = c(1/5, 1)));
#color.palette.div=rev(diverge_hcl(n=lev, h = c(260, 0), c = 80, l = c(30, 90), power = 1.5));
#filled.contour(f$x, f$y, z=f$z, nlevels=N,  col=color.palette.terrain);

filled.contour(f$x, f$y, z=f$z, nlevels=N,  col=color.palette.terrain, plot.axes={ points(exprs(filtdat$V1[,"FSC-A"]),exprs(filtdat$V1[,"SSC-A"]),col="red"); axis(1); axis(2)}, xlab="FSC-A", ylab="SSC-A",cex.lab=1.4)
##filled.contour(f$x, f$y, z=f$z, nlevels=N,  col=color.palette.terrain, plot.axes={ points(exprs(filtdat$V1[,"FSC-A"]),exprs(filtdat$V1[,"SSC-A"]),col=rgb(0,0,0,0)); axis(1); axis(2)}, xlab="FSC-A", ylab="SSC-A",cex.lab=1.4)
dev.off();

txtfilename=paste(txtoutname,fraccount,lengthfrac,sep="_");
txtfilename=paste(txtfilename,"txt",sep=".");
###write.table(exprs(filtdat$V1),txtfilename); 

meanexp=mean(exprs(filtdat$V1[,"FITC-A"]));
sdexp=sd(exprs(filtdat$V1[,"FITC-A"]));

print(lengthfrac);
print(meanexp);
print(sdexp);

lengthfrac=lengthfrac+0.25;

}

print("------------------------------------------------------------");

fraccount=fraccount+0.05;

}
sink();
