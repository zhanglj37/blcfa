
write_mplus_bayes<-function(varnames,usevar,myModel,filename,sigpsx_list)
{
if(file.exists("blcfa_bayes.inp"))
{
 file.remove("blcfa_bayes.inp")
}

SIGPSX=sigpsx_list$SIGPSX

## The variable name part is split into multiple lines in order to be smaller than the 90 characters constrained by Mplus.
if (length(varnames)%%10 == 0)
{
	var_row_num<- length(varnames)%/%10
}else{
	var_row_num<- length(varnames)%/%10+1
}

for (i in 1:var_row_num)
{
	if (i == var_row_num)
	{
		tempnum = length(varnames)-10*(i-1)
		if (tempnum == 1)
		{
			assign(paste0("combine_names",i),varnames[(i-1)*10+1])		
		}else if (tempnum ==2)
		{
			assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))		
		}else{
			assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))
			for (j in 3:tempnum)
			{
				assign(paste0("combine_names",i),paste(get(paste0("combine_names",i)),varnames[(i-1)*10+j]))
			}	
		}			
	}else{
		assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))
		for (j in 3:10)
		{
			assign(paste0("combine_names",i),paste(get(paste0("combine_names",i)),varnames[(i-1)*10+j]))
		}
	}
	
}

if (length(varnames)%%10 == 0)
{
	usevar_row_num<- length(usevar)%/%10
}else{
	usevar_row_num<- length(usevar)%/%10+1
}

for (i in 1:usevar_row_num)
{
	if (i == usevar_row_num)
	{
		tempnum = length(usevar)-10*(i-1)
		if (tempnum == 1)
		{
			assign(paste0("combine_usenames",i),usevar[(i-1)*10+1])		
		}else if (tempnum ==2)
		{
			assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))		
		}else{
			assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))
			for (j in 3:tempnum)
			{
				assign(paste0("combine_usenames",i),paste(get(paste0("combine_usenames",i)),usevar[(i-1)*10+j]))
			}	
		}			
	}else{
		assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))
		for (j in 3:10)
		{
			assign(paste0("combine_usenames",i),paste(get(paste0("combine_usenames",i)),usevar[(i-1)*10+j]))
		}
	}
	
}

## Model part: replace # with !, replace + with spaces, replace =~ with BY & spaces
model1<-gsub("#","!",myModel)
model2<-gsub("\\+","",model1)
model3<-gsub("\\=~"," BY ",model2)
model4<-gsub("\n",";\n\t",model3)


## generate input file
cat(
	"TITLE: Bayesian Lasso CFA\n", #_fix211_
	"DATA: FILE = ", filename, ";",
	"\n",
	"VARIABLE:\n",
	"NAMES = ",
	file = paste("blcfa_bayes.inp", sep = ''), append = T)

if (var_row_num == 1)
{  
	cat(get(paste0("combine_names",var_row_num)),";\n",
		file = paste("blcfa_bayes.inp", sep = ''), append = T)
	
}else{
	for (i in 1:(var_row_num-1))
	{
		cat(get(paste0("combine_names",i)),"\n\t",
			file = paste("blcfa_bayes.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_names",var_row_num)),";\n",
			file = paste("blcfa_bayes.inp", sep = ''), append = T)
}	
	
cat(
	"USEV = ",
	file = paste("blcfa_bayes.inp", sep = ''), append = T)

if (usevar_row_num == 1)
{
	cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
		file = paste("blcfa_bayes.inp", sep = ''), append = T)
}else{
	for (i in 1:(usevar_row_num-1))
	{
		cat(get(paste0("combine_usenames",i)),"\n\t",
			file = paste("blcfa_bayes.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
			file = paste("blcfa_bayes.inp", sep = ''), append = T)
} 
 
cat(
	"ANALYSIS:\n\t",
	"ESTIMATOR = BAYES;\n\t",
	"ALGORITHM=GIBBS(RW); \n\t",
	"PROC = 2;\n\t",
	"BITERATIONS = (10000);\n",
	"MODEL:\n\t",
	model4,
	"\n\n\t",
	file = paste("blcfa_bayes.inp", sep = ''), append = T)
 
if (SIGPSX[1] != 0)
{
	sigpsxname<-sigpsx_list$sigpsxname
	for (i in 1:length(sigpsxname))
	{
		cat(sigpsxname[i],";\n\t",
			file = paste("blcfa_bayes.inp", sep = ''), append = T)
	}
}

cat(
	"\n",
	"OUTPUT: TECH1  TECH8  STDY;\n",
	"PLOT: TYPE= PLOT2;\n",
	file = paste("blcfa_bayes.inp", sep = ''), append = T)
 
## run
runModels()

if(file.exists("Mplus Run Models.log"))
{
 file.remove("Mplus Run Models.log")
}

}