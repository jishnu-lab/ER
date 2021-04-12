#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

require(huge)
source('precisionMat2txt.R')

if (length(args) != 2) {
    stop("Missing arguments\nUsage is: Rscript runGlasso.R input_data output_file", call.=FALSE)
}

data <- read.delim(args[1])

out.glasso <- huge(as.matrix(data), method="glasso" , nlambda=30, lambda.min.ratio=0.1)

out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)

print(paste("optimal lambda:", out.select$opt.lambda, sep=" "))

precisionMat2txt(out.select$opt.icov, fname=args[2], nodes=colnames(data))

