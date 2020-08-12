get_chr_windows_hg19 <- function(wind){
  
  # number of bases in chromosomes
  chr_base<-list(
    c(249250621, 243199373, 198022430, 191154276,
      180915260, 171115067, 159138663, 146364022,
      141213431, 135534747, 135006516, 133851895,
      115169878, 107349540, 102531392, 90354753,
      81195210, 78077248, 59128983, 63025520,
      48129895, 51304566, 155270560, 59373566
      ),
    c(as.character(seq(1, 22, 1)), "X", "Y"))
  
  chr_windows <- data.frame(chr=0, from=0, to=0)
  chr_windows <- chr_windows[-1, ]
  
  # for each chromosome
  for (i in 1:length(chr_base[[2]])){
    
    # print(paste("Chromosome ",chr_base[[2]][i]))
    # generate windows boundaries including right border
    a <- seq(0, chr_base[[1]][i], wind)
    a <- append(a, chr_base[[1]][i])
    a <- unique(a)
    
    from <- a[1:length(a) - 1]
    to <- a[2:length(a)]
    
    new_row <- as.data.frame(cbind(rep(chr_base[[2]][i], length(a) - 1), from,to))
    names(new_row) <- names(chr_windows)
    chr_windows <- rbind(chr_windows, new_row)
    
  }
  chr_windows$from <- as.numeric(as.character(chr_windows$from))
  chr_windows$to <- as.numeric(as.character(chr_windows$to))
  return(chr_windows)
}


# set window size
window_size <- 10000
chr_windows <- get_chr_windows_hg19(window_size)

