#source("~/Work/AManal/MCMCglmm_2.07/inst/doc/Figures/buildCN.R")
alone=FALSE
lectures<-1:9
forCRAN=FALSE
JSS<-FALSE
LINUX=TRUE
options(width=80)

UPpath="~/Work/AManal/UP_course/Tex/"
MCpath="~/Work/AManal/MCMCglmm_2.07/inst/doc/"

library(MCMCglmm)
library(MCMCpack)
library(R2WinBUGS)
library(kinship)

data(BTdata)

if(alone==FALSE){
	setwd(MCpath)	
	for(i in lectures){
	   system(paste("rm Lecture", i, "*", sep=""))
	}
	system("rm CourseNotes.pdf")
	system("rm CourseNotes.tex")
	system("rm CourseNotes.log")
	system("rm CourseNotes.aux")
	system("rm CourseNotes.blg")
}  # remove files from master


for(i in lectures){
	
  UPpath_tmp<-paste(UPpath,  "Lecture", i, "/", sep="")

  setwd(UPpath_tmp)
  
  system("rm *.pdf")
  system("rm *.eps")
  system("rm *.aux")
  system("rm *.log")  # remove old files

  Rnw<-paste(UPpath_tmp, "Lecture", i, ".Rnw", sep="")
  Tex<-paste(UPpath_tmp, "Lecture", i, ".tex", sep="")
  Bib<-paste(UPpath_tmp, "Lecture", i, sep="")
  Pdf<-paste(UPpath_tmp, "Lecture", i, ".pdf", sep="")

  if(alone){
  	 system(paste("ex", Rnw, "-s +/alonefalse +:s/alonefalse/alonetrue +:x"))
  }else{ 
     system(paste("ex", Rnw, "-s +/alonetrue +:s/alonetrue/alonefalse +:x"))
  } 
  		
  Sweave(Rnw)
  
#	system(paste("ex",Tex, "+:s/alonefalse/alonetrue +:x"))

  if(alone){
  	if(LINUX){
	  system(paste("cp ", MCpath, "JarLib.bib ", UPpath_tmp, sep=""))

	  system(paste("pdflatex", Tex))
	  system(paste("bibtex", Bib))
	  system(paste("pdflatex", Tex))
	  system(paste("pdflatex", Tex))
	
	  if(length(lectures)==1){system(paste("evince", Pdf, "&"))}

       }else{

	  system(paste("cp ", MCpath, "JarLib.bib ", UPpath_tmp, sep=""))

	  system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	  system(paste("~/Library/TeXShop/bin/bibtexc", Bib))
	  system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	  system(paste("~/Library/TeXShop/bin/pdflatexc", Tex))
	
	  if(length(lectures)==1){system(paste("open -a Preview", Pdf))}


       }
    
  }else{

	system(paste("cp *.tex", MCpath))
	system(paste("cp *.pdf", MCpath))  # copy tex and pdf files over to master
	
  }	
  
}
			   
if(alone==FALSE){	
	setwd(MCpath)
	system(paste("cp ", MCpath, "Figures/CourseNotes.Rnw ", MCpath, sep="")) # copy master file back out of figures
	Sweave("CourseNotes.Rnw") 

    if(LINUX){
	system("pdflatex CourseNotes.tex")
	system("bibtex CourseNotes")
	system("pdflatex CourseNotes.tex")
	system("pdflatex CourseNotes.tex")
        system("evince CourseNotes.pdf&")
    }else{
	system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
	system("~/Library/TeXShop/bin/bibtexc CourseNotes")
	system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
	system("~/Library/TeXShop/bin/pdflatexc CourseNotes.tex")
        system("open -a Preview CourseNotes.pdf")
    }
}

if(forCRAN){
	system(paste("cp ", MCpath, "/CourseNotes.Rnw ", MCpath, "/Figures", sep=""))
	system(paste("cp ", MCpath, "/CourseNotes.tex ", MCpath, "/CourseNotes.Rnw", sep=""))
	system(paste("rm ", MCpath, "/CourseNotes.tex", sep=""))
#	system(paste("rm ", MCpath, "/JSS.bst", sep=""))
#	system(paste("rm ", MCpath, "/JSS.cls", sep=""))
}

if(JSS){
	setwd("~/Desktop/Work/Tex")
	Sweave("article.Rnw")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
	system("~/Library/TexShop/bin/bibtexc article")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
	system("~/Library/TexShop/bin/pdflatexc article.tex")
    system("open -a Preview article.pdf")
	if(redoOverview){
	setwd(MCpath)
	Sweave("Overview.Rnw")
	system("~/Library/TexShop/bin/pdflatexc Overview.tex")
	system("~/Library/TexShop/bin/bibtexc Overview")
	system("~/Library/TexShop/bin/pdflatexc Overview.tex")
	system("~/Library/TexShop/bin/pdflatexc Overview.tex")
    system("open -a Preview Overview.pdf")
	system(paste("cp ", MCpath, "/Overview.Rnw ", MCpath, "/Figures", sep=""))
	system(paste("cp ", MCpath, "/Overview.tex ", MCpath, "/CourseNotes.Rnw", sep=""))	
	}
}

