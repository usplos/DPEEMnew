
library(tidyverse)
library(rio)

seq_evermax = function(x){
  Log = logical(length = length(x))
  Log[1] = T
  if(length(x) > 1){
    for (ii in 2:length(x)) {
      Log[ii] = ifelse(x[ii] > max(x[1:(ii-1)]), T, F)
    }
  }
  return(Log)
}

seq_sacdir = function(x){
  Dir = character(length = length(x))
  if(length(x) == 1){
    Dir[1] = 'NO'
  }else{
    Dir = c(ifelse(x[-1] > x[-length(x)],'fore','back'),'NO')
  }
  return(Dir)
}

seq_locatetimes = function(x){
  times = numeric(length = length(x))
  target = unique(x[x>0])
  for(ii in target){
    xnew = which(x==ii)
    if(length(xnew)==1) {
      iitimes = 1
    }else{
      iitimes = c(1)
      for(jj in 2:(length(xnew))){
        if(xnew[jj] - xnew[jj-1] == 1){
          iitimes[jj] = iitimes[jj-1]
        }else{
          iitimes[jj] = iitimes[jj-1]+1
        }
      }
    }
    times[which(x==ii)] = iitimes
  }
  times
}

seq_regressionfrom = function(sacdir,roi){
  regression = rep('NO',times=length(sacdir))
  if(length(sacdir) > 1){
    for (rr in 2:length(sacdir)) {
      regression[rr] = ifelse(roi[rr]>0 ,
                              ifelse(roi[rr] == roi[rr-1],
                                     regression[rr-1],
                                     ifelse(sacdir[rr-1] == 'fore',
                                            'Left', 'Right')),
                              'NO')
    }
  }
  regression
}

funDA12csv <-function(workdir){
  # convert the DA1 files to csv files
  # the DA1 files should be names as Sub*.DA1
  # the DA1 files should be in current working directory
  setwd(workdir)
  DA1filename = dir(pattern = '[Ss]ub')
  DA1filename = DA1filename[grep('*.DA1',DA1filename)]
  if(length(DA1filename) == 0)
  {cat('Error: there is no DA1 file\n')}else{
    for (ida1 in DA1filename) {
      icsv = paste(substr(ida1,1,nchar(ida1)-4),'.csv',sep = '')
      file.copy(ida1, icsv)
      cat(ida1,'is done\n')
    }
    cat('DA1 to csv convert is done!\n\n')

  }
}

funpreprocess <-function(workdir, FDMax = 1000, FDMin = 80, ROIfilename, outputdir){
  # This function is aimed to tide the csv files of each subjects.
  # Files of all subjects will be tidied in to FTtotal.csv in which one line contains the basic information about one fixation point.
  # FTtotalA.csv will then be produced based on FTtotal.csv in which several variables is appended inclusing the final coordinate (xcoor + ycoor*160), whether it is the first-pass fixation, fixation duration, whether it is in the ROI, saccade direction from current fixation point.
  # FTtotalAS.csv will then be produced based on FTtotalA.csv in which fixations whose duration located within the limitation will be remained.
  # FTtotalASR.csv will then be produced based on FTtotalAS.csv in which only fixations locating in ROI are remained.
  # FTtotalASRpt.csv will then be produced based on FTtoalAS.csv in which another variable is appended which indicats how many times the current fixation point has passed the ROI if it is located in ROI.
  # FTtotalASRptReg.csv will finally be produced based on FTtotalASRptReg.csv in which another variable is appended which indicates whether the current saccade to ROIs is from right or left to ROIs.
  # workdir - a character string indicating where the csv files are
  # FDMax - numeric indicating the longest duration of fixation point you want to remain
  # FDMin - numeric indicating the shortest duration of fixation point you want to remain
  # ROIfilename - a character string indicating the name ROI location file
  # outputdir - a character string indicating where you want to put out the FTtotal*.csv files
  #workdir = getwd();FDMax = 1000;FDMin = 80;ROIfilename = 'ROI.xlsx';outputdir = getwd()
  setwd(workdir)
  cat('Now preprocessing is done...\n')

  csvfilename = dir(pattern = '[Ss]ub')
  csvfilename = csvfilename[grep('*.csv',csvfilename)]
  sub0 = c()
  item0 = c()
  cond0 = c()
  xcoor0 = c()
  ycoor0 = c()
  Tstart0 = c()
  Tend0 = c()

  for(i in csvfilename){
    tempsubfile = read.csv(i, sep = ' ', header = F)
    if(ncol(tempsubfile) <2){
      tempsubfile = read.csv(i, sep = '\t', header = F)
    }
    for(j in 1:nrow(tempsubfile)){
      templine = tempsubfile[j,]
      templine = templine[,!is.na(templine)]
      if(ncol(templine) > 8 | templine[[3]] > 199){
        if(templine[[3]] < 199){
          numberfixation = (ncol(templine)-8)/4

          for(k in 1:numberfixation){
            sub0 = c(sub0, substr(i,1,nchar(i)-4))
            cond0 = c(cond0, templine[[2]])
            item0 = c(item0, templine[[3]])
            xcoor0 = c(xcoor0, templine[[4*k+5]])
            ycoor0 = c(ycoor0, templine[[4*k+6]])
            Tstart0 = c(Tstart0, templine[[4*k+7]])
            Tend0 = c(Tend0, templine[[4*k+8]])
          }
        }else{
          numberfixation = ncol(templine)/4
          for(k in 1:numberfixation){
            sub0 = c(sub0, substr(i,1,nchar(i)-4))
            cond0 = c(cond0, templine[j-1,][[2]])
            item0 = c(item0, tempsubfile[j-1,][[3]])
            xcoor0 = c(xcoor0, templine[[4*k-3]])
            ycoor0 = c(ycoor0, templine[[4*k-2]])
            Tstart0 = c(Tstart0, templine[[4*k-1]])
            Tend0 = c(Tend0, templine[[4*k]])
          }
        }
      }
    }
    cat("\015")
    cat(i,' has been done!!',"\n")
  }

  FTtotal = tibble(sub0,cond0,item0,xcoor0,ycoor0,Tstart0,Tend0)
  naposition = which(is.na(FTtotal$cond0))
  for(i in naposition){
    FTtotal$cond0[i] = FTtotal$cond0[i-1]
  }
  largeitemp = which(FTtotal$item0 >= 199)
  for(i in largeitemp){
    FTtotal$item0[i] = FTtotal$item0[i-1]
  }

  write.csv(FTtotal,paste(outputdir,'FTtotal.csv',sep = '/'),row.names = F, quote = F)
  #FTtotal = read.csv(paste(outputdir,'FTtotal.csv',sep = '/'), stringsAsFactors = F)
  cat('FTtotal.csv has been produced','\n')

  FTtotal['finalcoor'] = FTtotal$xcoor0+FTtotal$ycoor0*160 # finalcoor

  ROI = import(ROIfilename)
  FTtotal['ROI0'] = 0
  for (ii in unique(ROI$item)) {
    ROIii = ROI %>% filter(item == ii)
    for(cc in unique(ROIii$condition)){
      ROIiicc = ROIii %>% filter(condition == cc) %>% arrange(roistart)
      for(rr in 1:nrow(ROIiicc)){
        FTtotal$ROI0[FTtotal$cond0 == cc &
                       FTtotal$item0 == ii &
                       FTtotal$finalcoor %in%
                       c(ROIiicc$roistart[rr] : (ROIiicc$roistart[rr] + ROIiicc$len[rr] - 1))] = ROIiicc$label[rr]
      }
    }


  } # ROI

  FTtotal['FFT'] = FTtotal$Tend0 - FTtotal$Tstart0 # FFT

  FTtotal = FTtotal %>% group_by(sub0, cond0, item0) %>% mutate(ffd0 = seq_evermax(finalcoor)) # ffd0

  FTtotal = FTtotal %>% group_by(sub0, cond0, item0) %>% mutate(sacdir = seq_sacdir(finalcoor)) # sacdir


  FTtotalA = FTtotal
  write.csv(FTtotalA, paste(outputdir, 'FTtotalA.csv',sep = '/'), row.names = F, quote = F)
  cat('FTtotalA.csv has been produced','\n')


  FTtotalAS = FTtotalA[FTtotalA$FFT %in% FDMax:FDMin,]
  write.csv(FTtotalAS, paste(outputdir, 'FTtotalAS.csv', sep = '/'),row.names = F, quote = F)
  cat('FTtotalAS.csv has been produced','\n')


  FTtotalASR = FTtotalAS[FTtotalAS$ROI0 > 0,]
  write.csv(FTtotalASR, paste(outputdir,'FTtotalASR.csv',sep = '/'), row.names = F, quote = F)
  cat('FTtotalASR.csv has been produced','\n')
  cat('Preprocessing has been done','\n\n')


  cat('Calculating the ROI passing times...\n')

  #FTtotalAS = read.csv(paste(outputdir, 'FTtotalAS.csv', sep = '/'), stringsAsFactors = F)
  FTtotalASRpt = FTtotalAS %>% group_by(sub0, item0, cond0) %>%
    mutate(passtimes = seq_locatetimes(ROI0))
  write.csv(FTtotalASRpt, paste(outputdir, 'FTtotalASRpt.csv', sep = '/'), quote = F, row.names = F)
  cat('FTtotalASRpt.csv has been done','\n\n')

  cat('Calculating whether regression was existed...\n')
  #FTtotalASRpt = read.csv(paste(outputdir,'FTtotalASRpt.csv', sep = '/'), stringsAsFactors = F)
  FTtotalASRptReg = FTtotalASRpt %>% group_by(sub0, item0, cond0) %>%
    mutate(regressionfrom = seq_regressionfrom(sacdir, ROI0))
  write.csv(FTtotalASRptReg, paste(outputdir, 'FTtotalASRptReg.csv', sep = '/'), row.names = F, quote = F)
  cat('FTtotalASRptReg.csv has been done','\n\n')

}

## temporal measurements
funROITTFFD <-function(workdir,outputdir){
  # calculate the first fixation duration and total time in ROI on each trial in each subject
  # based on FTtoalASR.csv
  # outputdir - a character string indicating where you want to put out the result file
  # FTtotalASR.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the first fixation duration and total time...\n')

  read_csv('FTtotalASRptReg.csv') %>%
    mutate(Sub = sub0, Cond = cond0, Item = item0,ROI = ROI0) %>%
    filter(ROI > 0) %>%
    group_by(Sub, Cond, Item, ROI) %>%
    summarise(TotalTime = sum(FFT)) %>%
    write_csv(paste0(outputdir,'/ROITT.csv'))

  read_csv('FTtotalASRptReg.csv') %>%
    mutate(Sub = sub0, Cond = cond0, Item = item0,ROI = ROI0) %>%
    filter(ROI0 > 0, passtimes == 1) %>%
    group_by(Sub, Cond, Item, ROI) %>%
    mutate(Count = 1:length(FFT)) %>%
    filter(Count == 1) %>%
    rename(FFD = FFT) %>%
    select(Sub, Cond, Item,ROI, FFD) %>%
    write_csv(paste0(outputdir,'/ROIFFD.csv'))


  cat('Totaltime and first fixation duration have been done','\n\n')

}

funROIgazeduration <-function(workdir,outputdir){
  # calcualte the gaze duration on ROI
  # based on FTtotalASRptReg.csv file
  # outputdir - a character string indicating where you want to put out the result files
  # FTtotalASRptReg.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the first pass time...\n')

  read_csv('FTtotalASRptReg.csv') %>%
    mutate(Sub = sub0, Cond = cond0, Item = item0, ROI = ROI0) %>%
    filter(passtimes == 1) %>%
    group_by(Sub, Cond, Item, ROI) %>%
    summarise(GazeDuration = sum(FFT)) %>%
    write_csv(paste0(outputdir,'/ROIGazeDuration.csv'))

  cat('ROI gaze duration has been done','\n\n')

}

funROIsecondFT <-function(workdir,outputdir){
  # calculate the duration of the second, third and four times passing the ROI
  # based on FTtotalASRepReg.csv
  # outputdir - a character string indicating where you want to put out the result file
  # FTtotalASRptReg.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the second pass time...\n')
  FTtotalASRptReg  = read_csv('FTtotalASRptReg.csv')

  #passtimes = FTtotalASRptReg %>% filter(sub0 == 'Sub_002', item0 == 3, ROI0==3) %>% .$passtimes
  #FFT = FTtotalASRptReg %>% filter(sub0 == 'Sub_002', item0 == 3, ROI0==3) %>% .$FFT
  secondpasstime = function(passtimes, FFT){
    ifelse(length(which(passtimes == 2)) == 0,0,
           FFT[which(passtimes == 2)] %>% sum())
  }
  thirdpasstime = function(passtimes, FFT){
    ifelse(length(which(passtimes == 3)) == 0,0,
           FFT[which(passtimes == 3)] %>% sum())
  }
  fourthpasstime = function(passtimes, FFT){
    ifelse(length(which(passtimes == 4)) == 0,0,
           FFT[which(passtimes == 4)] %>% sum())
  }

  FTtotalASRptReg %>% rename(Sub = sub0, Item = item0, Cond = cond0, ROI = ROI0) %>%
    filter(ROI > 0) %>%
    group_by(Sub, Item, Cond, ROI) %>%
    summarise(Secondpasstime = secondpasstime(passtimes, FFT),
              Thirdpasstime = thirdpasstime(passtimes, FFT),
              Fourthpasstime = fourthpasstime(passtimes, FFT)) %>%
    arrange(Sub, Item, Cond,ROI) %>%
    write.csv(paste(outputdir,'ROIsecondFT.csv', sep = '/'), row.names = F, quote = F)

  cat('ROI second fixation duration has been done','\n\n')

}

# integer measurements
funFTROInum <-function(workdir,outputdir){
  # calculate the fixation numbers and its proportion in ROI for each trial on each prticipant
  # based on FTtotalAS.csv
  # outputdir - a character string indicating where you want to put out the result file
  # FTtotalAS.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the number of fixation in ROI...\n')
  ROI = import('ROI.xlsx')
  ROIfixationnumber = function(df){
    # df = FTtotalASRptReg %>% filter(Sub == 'Sub_001', Item == 1)
    numberroi = ROI %>% filter(item == df$Item[1],
                               condition == df$Cond[1]) %>% nrow()
    Fixationnumber = c()
    Fixationproportion = c()
    for(rr in 1:numberroi){
      dfrr = df %>% filter(ROI == rr)
      Fixationnumber[rr] = nrow(dfrr)
      Fixationproportion[rr] = nrow(dfrr)/nrow(df)
    }

    tibble(Sub = rep(df$sub0[1], numberroi),
           Item = rep(df$item0[1], numberroi),
           Cond = rep(df$cond0[1], numberroi),
           ROI = 1:numberroi,
           Fixationnumber,
           Fixationproportion)
  }

  FTtotalASRptReg = read_csv('FTtotalASRptReg.csv') %>%
    mutate(Sub = sub0, Cond = cond0, Item = item0, ROI=ROI0)
  ROIfixationnumbertb = tibble()
  for(ss in unique(FTtotalASRptReg$Sub)){
    dfsub = FTtotalASRptReg %>% filter(Sub %in% ss)
    for(ii in unique(dfsub %>% filter(ROI > 0) %>% .$Item)){
      dfsubitem = dfsub %>% filter(Item == ii)
      for(cc in unique(dfsubitem$Cond)){
        dfsubitemcond = dfsubitem %>% filter(Cond==cc)
        ROIfixationnumbertb = bind_rows(ROIfixationnumbertb, ROIfixationnumber(df = dfsubitemcond))
      }
    }
  }

  ROIfixationnumbertb %>% arrange(Sub, Item, Cond, ROI) %>%
    write_csv(paste0(outputdir,'/ROIFTnum.csv'))
  cat('FTROInum.csv has been produced','\n\n')
}

funROIsaccadelength <-function(workdir,outputdir){
  # calculate the saccade length saccading into and out-of ROI when first time passing based on the direction of saccade from
  # based on FTtotalASRptReg.csv file
  # outputdir - a character string indicating where you want to put out the result file
  # FTtotalASRptReg.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the saccade length...\n')
  FTtotalASRptReg  = read_csv('FTtotalASRptReg.csv')

  saccadelength = function(df){
    lengthIntoROIleft = c()
    lengthIntoROIright = c()
    lengthOutfromROIleft = c()
    lengthOutfromROIright = c()
    for(rr in unique(df$ROI0)[which(unique(df$ROI0)>0)]){
      df2 = df %>% mutate(ROI0 = ifelse(ROI0 == rr & passtimes == 1, rr, 0))

      lengthIntoROIleft = c(lengthIntoROIleft,ifelse(which(df2$ROI0 == rr)[1] == 1, 0,
                                                     ifelse(df2$sacdir[which(df2$ROI0 == rr)[1] -1] == 'fore',
                                                            abs(df2$finalcoor[(which(df2$ROI0 == rr)[1])-1] - df2$finalcoor[which(df2$ROI0 == rr)[1]]),0)))
      lengthIntoROIright = c(lengthIntoROIright,ifelse(which(df2$ROI0 == rr)[1] == 1, 0,
                                                       ifelse(df2$sacdir[which(df2$ROI0 == rr)[1] -1] == 'back',
                                                              abs(df2$finalcoor[(which(df2$ROI0 == rr)[1])-1] - df2$finalcoor[which(df2$ROI0 == rr)[1]]),0)))

      lengthOutfromROIleft = c(lengthOutfromROIleft,ifelse(which(df2$ROI0 == rr) %>% .[length(.)] == nrow(df2),0,
                                                           ifelse(df2$sacdir[which(df2$ROI0 == rr) %>% .[length(.)]] == 'fore',
                                                                  abs(df2$finalcoor[which(df2$ROI0 == rr) %>% .[length(.)]] -
                                                                        df2$finalcoor[(which(df2$ROI0 == rr) %>% .[length(.)]) + 1]),0)))
      lengthOutfromROIright = c(lengthOutfromROIright,ifelse(which(df2$ROI0 == rr) %>% .[length(.)] == nrow(df2),0,
                                                             ifelse(df2$sacdir[which(df2$ROI0 == rr) %>% .[length(.)]] == 'back',
                                                                    abs(df2$finalcoor[which(df2$ROI0 == rr) %>% .[length(.)]] -
                                                                          df2$finalcoor[(which(df2$ROI0 == rr) %>% .[length(.)]) + 1]),0)))
    }

    tibble(Sub = rep(df$sub0[1], length(unique(df$ROI0)[which(unique(df$ROI0)>0)])),
           Item = rep(df$item0[1],length(unique(df$ROI0)[which(unique(df$ROI0)>0)])),
           Cond = rep(df$cond0[1], length(unique(df$ROI0)[which(unique(df$ROI0)>0)])),
           ROI = unique(df$ROI0)[which(unique(df$ROI0)>0)],
           lengthIntoROIleft,
           lengthIntoROIright,
           lengthOutfromROIleft,
           lengthOutfromROIright)
  }

  ROIsaccadelength = tibble()
  for(ss in unique(FTtotalASRptReg$sub0)){
    dfsub = FTtotalASRptReg %>% filter(sub0 %in% ss)
    for(ii in unique(dfsub %>% filter(ROI0 > 0) %>% .$item0)){
      dfsubitem = dfsub %>% filter(item0 == ii)
      for(cc in unique(dfsubitem$cond0)){
        dfsubitemcond = dfsubitem %>% filter(cond0==cc)
        ROIsaccadelength = bind_rows(ROIsaccadelength, saccadelength(df = dfsubitemcond))
      }
    }
  }

  ROIsaccadelength %>% arrange(Sub, Item, Cond, ROI) %>%
    write.csv(paste(outputdir,'ROIsaccadelength.csv', sep = '/'), row.names = F, quote = F,na = '')
  cat('ROI saccade length is done\n\n')

}

# logical mearsurements
funROIfixationprop <-function(workdir,outputdir){
  # calculate whether the ROI was focused when subject passed it first time
  # based on FTtotalASRptReg.csv file
  # outputdir - a character string where you want to put out the result file
  # FTtotalASRptReg.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the fixation proportion...\n')
  ROI = import('ROI.xlsx')
  firstpass = function(df){
    numberroi = ROI %>% filter(item %in% df$Item[1],
                               condition %in% df$Cond[1]) %>% nrow()
    fixationprop = c()
    for(ii in 1:numberroi){
      fixationprop[ii] = ifelse(length(df$ROI[which(df$ROI==ii)]) == 0,
                                0,
                                ifelse(isTRUE(df$ffd0[which(df$ROI==ii)][1]),
                                       1,0))
    }
    tibble(Sub = df$sub0[1] %>% rep(length(fixationprop)),
           Item = rep(df$item0[1], times=length(fixationprop)),
           Cond = rep(df$cond0[1], times=length(fixationprop)),
           ROI = 1:numberroi,
           FirstpassFixated = fixationprop,
           FirstpassSkipped = 1 - fixationprop)
  }

  FTtotalASRptReg = read_csv('FTtotalASRptReg.csv') %>%
    mutate(Sub = sub0, Cond = cond0, Item = item0, ROI = ROI0) %>%
    filter(ROI0 > 0)

  ROIfixationprop = tibble()
  for(ss in unique(FTtotalASRptReg$Sub)){
    dfsub = FTtotalASRptReg %>% filter(Sub %in% ss)
    for(ii in unique(dfsub$Item)){
      dfsubitem = dfsub %>% filter(Item == ii)
      for(cc in unique(dfsubitem$Cond)){
        dfsubitemcond = dfsubitem %>% filter(Cond==cc)
        ROIfixationprop = bind_rows(ROIfixationprop,
                                    firstpass(df = dfsubitemcond))
      }
    }
  }
  ROIfixationprop %>% arrange(Sub, Item, Cond, ROI) %>%
    write_csv(paste0(outputdir, '/ROIFixationProp.csv'))




  cat('ROI fixation proportion has been done','\n\n')

}

funROIregressionIn_Out <-function(workdir,outputdir){
  # calculate whether the ROI received regression and whether the first time landing in ROI was regression
  # based on FTtotalASRepReg.csv
  # outputdir - a character string indicating where you want put out the result file
  # FTtotalASRptReg.csv should be in current working directory
  setwd(workdir)
  cat('Calculating the ROI regression ...\n')

  ROI = import('ROI.xlsx')
  regressionin_out = function(df){
    numberroi = ROI %>% filter(item == df$item0[1],
                               condition == df$cond0[1]) %>% nrow()
    RegressionIn = c()
    RegressionInfirstpass = c()
    RegressionOut = c()
    RegressionOutfirstpass = c()
    for(rr in 1:numberroi){
      RegressionIn[rr] = ifelse(nrow(df %>% filter(ROI0 == rr)) == 0, 0,
                                ifelse(length(which(df %>% filter(ROI0 == rr) %>% .$regressionfrom == 'Right')) > 0,
                                       1,0))
      RegressionInfirstpass[rr] = ifelse(nrow(df %>% filter(ROI0==rr, passtimes == 1)) == 0, 0,
                                         ifelse(length(which(df %>% filter(ROI0==rr, passtimes == 1) %>% .$regressionfrom == 'Right')) > 0,
                                                1,0))
      RegressionOut[rr] = ifelse(nrow(df %>% filter(ROI0 == rr)) == 0, 0,
                                 ifelse(length(which(df %>% filter(ROI0 == rr) %>% .$sacdir == 'back')) > 0,
                                        1,0))
      RegressionOutfirstpass[rr] = ifelse(nrow(df %>% filter(ROI0==rr, passtimes == 1)) == 0, 0,
                                          ifelse(length(which(df %>% filter(ROI0==rr, passtimes == 1) %>% .$sacdir == 'back')) > 0,
                                                 1,0))
    }

    tibble(Sub = rep(df$sub0[1], numberroi),
           Item = rep(df$item0[1], numberroi),
           Cond = rep(df$cond0[1], numberroi),
           ROI = 1:numberroi,
           RegressionIn,
           RegressionInfirstpass,
           RegressionOut,
           RegressionOutfirstpass)
  }

  FTtotalASRptReg = read_csv('FTtotalASRptReg.csv')
  ROIregression = tibble()
  for(ss in unique(FTtotalASRptReg$sub0)){
    dfsub = FTtotalASRptReg %>% filter(sub0 %in% ss)
    for(ii in unique(dfsub %>% filter(ROI0 > 0) %>% .$item0)){
      dfsubitem = dfsub %>% filter(item0 == ii)
      for(cc in unique(dfsubitem$cond0)){
        dfsubitemcond = dfsubitem %>% filter(cond0==cc)
        ROIregression = bind_rows(ROIregression, regressionin_out(df = dfsubitemcond))
      }
    }
  }
  ROIregression %>% arrange(Sub, Item, Cond, ROI) %>%
    write.csv(paste(outputdir,'ROIregressionIn_Out.csv', sep = '/'),
              row.names = F, quote = F)

  cat('ROI regression In and Out has been done','\n\n')

}

# integrate process
funIntegrate <-function(workdir = getwd(),
                        outputdir = workdir,
                        FDMax, FDMin,ROIfilename1,

                        DA1_to_csv = T,preprocess = T,

                        TTFFD = T,GazeDuration = T,SecondPassTime = T,

                        FTnum = T,SaccadeLength = T,

                        Regression = T,FixationProportion = T,

                        DataIntegrate = T){
  # integrate the DA1 converting to csv, preprocessing and measures extracting
  # rio package is need, but dont worry about it because this function will check whether this package has been downloaded and will download it if not.
  # workdir - a character string indicating where the files you wang to use first are, defaultly the current working directory
  # outputdir - a character string indicating where you want to put out all the preprocessing and result files, can be a vector if you have more than one ROI
  # FDMax - numeric indicating the longest duration of fixation point you want to remain
  # FDMin - numeric indicating the shortest duration of fixation point you want to remain
  # ROIfilename1 - a character string indicating the name of the ROI location file (can be a vector if you have more than one ROI)
  # DA1_to_csv - logical indicating whether to convert the DA1 files
  # preprocess - logical indicating whether to preprocess
  # TTFFD - logical indicating whether to calculating total time and first fixation duration in ROI
  # FTnum - logical indicating whether to calculating the fixation numbers and its proportion in ROI
  # GazeDuration - logical indicating whether to calculating the gaze duration through ROI
  # Regression - logical indicating whether to calculating data regarding regression in ROI
  # SaccadeLength - logical indicating whether to calculating the saccade length to and from ROI
  # SecondPassTime - logical indicating whether to calculating the second pass time through ROI
  # FixationProportion - logical indicating whether to calculating the fixation proportion when passing ROI first time
  # DataIntegrate - logical whether to integrate all the result files into one result file

  # check whether the rio package has been downloaded
  if(!require(rio)) install.packages('rio')
  if(!require(tidyverse)) install.packages('tidyverse')
  library(rio);library(tidyverse)

  if(DA1_to_csv == T) funDA12csv(workdir = workdir)

  if(preprocess == T) funpreprocess(workdir = workdir,
                                    FDMax = FDMax, FDMin = FDMin,
                                    ROIfilename = ROIfilename1,
                                    outputdir = outputdir)

  checkpreprocess = dir(pattern = 'FTtotalASRptReg.csv')
  if(length(checkpreprocess) == 1){

    if(TTFFD == T) funROITTFFD(workdir = workdir, outputdir = outputdir)
    if(GazeDuration == T) funROIgazeduration(workdir = workdir, outputdir = outputdir)
    if(SecondPassTime == T) funROIsecondFT(workdir = workdir, outputdir = outputdir)

    if(FTnum == T) funFTROInum(workdir = workdir, outputdir = outputdir)
    if(SaccadeLength == T) funROIsaccadelength(workdir = workdir, outputdir = outputdir)

    if(Regression == T) funROIregressionIn_Out(workdir = workdir, outputdir = outputdir)
    if(FixationProportion == T) funROIfixationprop(workdir = workdir, outputdir = outputdir)

    if(DataIntegrate == T){
      setwd(workdir)
      datafilename = dir(pattern = 'ROI[a-zA-Z]')
      datafilename = datafilename[grep('.csv', datafilename)]
      datafilename = datafilename[which(datafilename != "ROITotal.csv")]

      if(length(datafilename) == 0){
        cat('Error: there is no data file!')
      }
      if(length(datafilename) == 1){
        cat('There is only one file, so there is no need to integrate')
      }
      if(length(datafilename) > 1){
        ROItotal = import(datafilename[1])
        for (ii in 2:length(datafilename)) {
          df2 = import(datafilename[ii])
          ROItotal = full_join(ROItotal,df2, by = c('Sub','Item','Cond','ROI'))
        }

        for(cc in 5:ncol(ROItotal)){
          ROItotal[[cc]][which(is.na(ROItotal[[cc]]))] = 0
        }
        export(ROItotal, paste(outputdir, 'ROITotal.csv', sep = '/'))
        cat('Data Integrate is done!\n\n')

      }

      cat('Congratulations!!!\n\n')
    }
  }else{
    cat('\nError:You have not done the preprocess analysis, please be sure that you have done that!\n\n')
  }
}

funGUI <-
  function()
  {
    
    if(sum(unique(installed.packages()[,c('Package')] %in% 'tools')) == 0)
    {install.packages('tools')}
    if(sum(unique(installed.packages()[,c('Package')] %in% 'fgui')) == 0)
    {install.packages('fgui')}
    library(fgui)
    
    res = gui(funIntegrate,
              argOption = list(DA1_to_csv = c('T','F'), preprocess = c('T','F'), TTFFD = c('T','F'),
                               FTnum = c('T','F'), GazeDuration = c('T','F'),
                               Regression = c('T','F'), SaccadeLength = c('T','F'),
                               SecondPassTime = c('T','F'), FixationProportion = c('T','F'),
                               SkipRate = c('T', 'F'),
                               DataIntegrate = c('T','F')),
              #argEdit = list(FDMax = NULL, FDMin = NULL),
              title = 'advanced DPEEM')
    
  }
