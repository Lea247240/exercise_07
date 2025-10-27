rm(list=ls())

#Cesta:
setwd('V:/MPA_PRG/exercise_07')
setwd('D:/VUT/4-5rocnik/moje/MPA-PRG/exercise_07') # home

# ..........................................................................................................
# Motif Search

## The Brute-Force Motif Search
### Task 1
# In R, create a function `Score()`, that calculates the score for a consensus string.

# Input:
  #An array of starting indexes.
  #`DNAStringSet` object of sequences (for example file `seq_score.fasta`).
  #Motif length.

# Output:
  # The score for the consensus string.

# call libraries
library(Biostrings)

# read seq_score.fasta
dna <- readDNAStringSet('seq_score.fasta', format = 'fasta')
dna

#initialization input
s <- c(1,30,8,19,44)
i <- 6

#....................................
# pokus
length(s)
s[1]+i

dna[[1]][1:3] # first seq 1-3 letters

dna[[1]][s[1]:(s[1]+i)]
dna[[1]][8:14]

# vybrat usek 
matrix <- c()

for (j in 1:length(s)){
  matrix <- c(matrix, dna[[j]][s[1]:(s[1]+i)])
}

matrix
matrix[[1]][1]

#not necessary
# vypocet cetnosti nukleotidu v sekvenci: vsechny <- alphabetFrequency(d) 
#                                       : pocetA <- letterFrequency(d, 'A')
# seznam_sek <- DNAStringSet(c('GA','GAN', 'NGN' ))
#pocetA <- letterFrequency(dna[1], 'A')



#vybrane useky podle zacatku indexu a dane do vektoru ............................... without function(debugging)
alignment_matrix <- c()
  
for (j in 1:length(s)){
  print(dna[[j]])
  for (k in 1:length(i))
    print(dna[[j]][s[j]:(s[j]+(i-1))])
    alignment_matrix <- c(alignment_matrix, dna[[j]][s[1]:(s[1]+(i-1))])
  
}
alignment_matrix

# prevedeni na stringset
alignment_matrix <- DNAStringSet(alignment_matrix)
alignment_matrix

# frequences profile
frequency_profile <- consensusMatrix(alignment_matrix)

#consensus string
consensus_string <- consensusString(frequency_profile)
consensus_string

#consensus score
score <- 0
for (j in 1:i){
  score_sloupec <- max(frequency_profile[,j])
  print(score_sloupec)
  score <- score + score_sloupec
}

score

#  ...................................................................................... all function (work)
Score <- function(s,dna,i){
  
  # create matrix alignment (find part seq tech motifs for start index (s) a get this together)
  alignment_matrix <- c()
  
  for (j in 1:length(s)){
    print(dna[[j]])
    for (k in 1:length(i))
      print(dna[[j]][s[j]:(s[j]+(i-1))])
    alignment_matrix <- c(alignment_matrix, dna[[j]][s[1]:(s[1]+(i-1))]) }
    
  alignment_matrix <- DNAStringSet(alignment_matrix)
  
  # frequency profile
  frequency_profile <- consensusMatrix(alignment_matrix)
  
  #consensus string
  consensus_string <- consensusString(frequency_profile)
  
  #consensus score
  score <- 0
  for (j in 1:i){
    score_sloupec <- max(frequency_profile[,j])
    print(score_sloupec)
    score <- score + score_sloupec
  }
  
  return(score) 
}

library(Biostrings)
Score(s = c(1,30,8,19,44) ,dna = readDNAStringSet('seq_score.fasta', format = 'fasta'),i = 6)


### Task 2 # .......................................................................................................
# In R, create function `NextLeaf()` according to the following pseudocode.

# Input:
  # `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
  # `t` Number of sequences.
  # `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

# Output:
  # `s` An array of starting indexes that corresponds to the next leaf in the tree.

################################
#NextLeaf(s, t, k)
#   for i ← t to 1
#     if s[i] < k
#       s[i] ← s[i] + 1
#       return s
#     s[i] ← 1
#   return s
################################

# read seq_score.fasta
dna <- readDNAStringSet('seq_score.fasta', format = 'fasta')
dna

#initialization input
s <- c(1,30,8,19,44)
t <- 5
i <- 6
k <- ((length(dna[[1]]) - i) +1)


# all function 
NextLeaf <- function(s,t,k){
  for (i in t:1){
    if (s[i] < k){
      s[i] <- s[i] + 1
      return(s)}
    s[i] <- 1
    return(s)
  }
}

dna <- readDNAStringSet('seq_score.fasta', format = 'fasta')
NextLeaf(s=c(1,30,8,19,44),t=length(dna), k=((length(dna[[1]]) - i) +1))


#........................................................................................................ (check how it works)
### Task 3
# In R, create a function `BFMotifSearch()` according to the following pseudocode.

# Input:
  # `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
  # `t` Number of sequences.
  # `n` Length of each sequence.
  # `l` Motif length.

# Output:
  # `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

#############################################
#BFMotifSearch(DNA, t, n, l)
#1   s ← (1, 1, ... , 1)                    # !!!!!
#2   bestScore ← Score(s, DNA, l)
#3   while forever
#4     s ← NextLeaf(s, t, n − l + 1)
#5       if Score(s, DNA, l) > bestScore
#6         bestScore ← Score(s, DNA, l)
#7         bestMotif ← (s1, s2, . . . , st) # !!!!
#8       if s = (1, 1, . . . , 1)
#9         return bestMotif
############################################

# read seq_motif.fasta
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
DNA

#initialization input
s <- c(1,3,8)
t <- length(DNA)
l <- 3
n <- length(DNA[[1]])
k <- ((length(DNA[[1]]) - l) +1)

# all function 

BFMotifSearch <- function(DNA, t, n, l){
  s <- c(1,1,1)  
  bestScore <- Score(s, DNA, l)
  while (TRUE) {
    s <- NextLeaf(s, t, ((n-l)+1))
    if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- c(s1, s2, st)
    }
    if (s == c(1,1,1)){
      return(bestMotif)
    }
  }
}

BFMotifSearch(DNA = readDNAStringSet('seq_motif.fasta', format = 'fasta'), t = length(DNA) , n = length(DNA[[1]]), l=3)



## The Branch-and-Bound Motif Search .................................................................................
### Task 4
# In R, create a function `NextVertex()` according to the following pseudocode.

# Input:
  # `s` An array of starting indexes `s = (s1 s2 … st)`, where *t* is the number of sequences.
  # `i` Level of vertex.
  # `t` Number of sequences.
  # `k` `k = n - l + 1`, where `n` is length of sequences and `l` is motif length.

# Output:
  # `s` The next vertex in the tree.
  # Current level of vertex.

#####################################
#NextVertex(s, i, t, k)
#1   if i < t
#2     s[i + 1] ← 1
#3     return (s, i + 1)
#4   else
#5      for j ← t to 1
#6       if s[j] < k
#7         s[j] ← s[j] + 1
#8         return (s, j)
#9   return (s, 0)
###################################

# read seq_motif.fasta
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
DNA

#initialization input
s <- c(1,3,8)
t <- length(DNA)
l <- 3
n <- length(DNA[[1]])
k <- ((length(DNA[[1]]) - l) +1)

i <- 0 # level of vertex????what that mean

# all function

NextVertex <- function(s, i, t, k){
  if (i < t){
    s[i+1] <- 1
    return(s,i+1)
  }
  else {
    for (j in t:1){
      if (s[j] < k){
        s[j] <- s[j] + 1
        return(s,j)
      }
    }
  }
  return(s,0)
}

NextVertex(s=c(1,3,8), i=0, t=length(DNA),k=((length(DNA[[1]]) - l) +1))


### Task 5 ................................................................................................. (check at home)
# In R, create a function `ByPass()` according to the following pseudocode.

# Input:
  # `s = (s1 s2 … st)`; an array of starting indexes, where *t* is the number of sequences
  # `i`; level of vertex
  # `t`; number of DNA sequences
  # `k = n - l + 1`, where `n` is length of DNA sequences and `l` is motif length

# Output:
  # the next leaf after a skip of a subtree
  # current level of vertex

####################################
#ByPass(s, i, t, k)
#1   for j ← i to 1
#2     if s[j] < k
#3       s[j] ← s[j] + 1
#4       return (s, j)
#5   return (s, 0)
###################################

# read seq_motif.fasta
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
DNA

#initialization input
s <- c(1,3,8)
t <- length(DNA)
l <- 3
n <- length(DNA[[1]])
k <- ((length(DNA[[1]]) - l) +1)

i <- 0 # level of vertex????what that mean

# all function
ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(s,j)
    }
  }
  return(s,0)
}

ByPass(s=c(1,3,8), i=0, t=length(DNA),k=((length(DNA[[1]]) - l) +1))


### Task 6 ....................................................................................................
# In R, create a function `BBMotifSearch()` according to the following pseudocode.

# Input:
  # `DNA` `DNAStringSet` object of sequences (for example file `seq_motif.fasta`).
  # `t` Number of sequences.
  # `n` Length of each sequence.
  # `l` Motif length.
# Output:
  # `bestMotif` An array of starting positions for each sequence with the best score for the consensus string.

# Modify function `Score()` to calculate score for the consensus string of the first `i` sequences of `DNA`.

###############################################################################
#BBMotifSearch(DNA, t, n, l)
#1   s ← (1, ... , 1)
#2   bestScore ← 0
#3   i ← 1
#4   while i > 0
#5     if i < t
#6       optimisticScore ← Score(s, i, DNA, l) + (t - i) * l
#7       if optimisticScore < bestScore
#8         (s, i) ← ByPass(s, i, t, n - l + 1)
#9       else
  #10        (s, i) ← NextVertex(s, i, t, n − l + 1)
#11    else
 # 12      if Score(s, t, DNA, l) > bestScore
#13        bestScore ← Score(s, t, DNA, l)
#14        bestMotif ← (s1, s2, ... , st)
#15      (s, i) ← NextVertex(s, i, t, n − l + 1)
#16  return bestMotif
###############################################################################

# read seq_motif.fasta
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
DNA

#initialization input
s <- c(1,3,8)
t <- length(DNA)
l <- 3
n <- length(DNA[[1]])
k <- ((length(DNA[[1]]) - l) +1)
i <- 0

# all function

BBMotifSearch <- function(DNA, t, n, l){
  s <- c(1,1,1) # vektor 1 podle delky motivu l
  bestScore <- 0
  i <- 1
  while (i > 0){
    if (i<t){
      optimisticScore <- Score(s, i, DNA, l) + ((t - i) * l)
      if (optimisticScore < bestScore){
        s,i <- ByPass(s, i, t, n - l + 1) # !!!!!! chyba
      }
      else {
        s, i <- NextVertex(s, i, t, ((n-l) + 1)) # !!!!!! chyba
      }
    }
    else {
      if (Score(s, t, DNA, l) > bestScore){
        bestScore <- Score(s, t, DNA, l)
        bestMotif <- (s1, s2, ... , st) # !!!!!! chyba
      }
      s, i <- NextVertex(s, i, t, ((n-l) + 1))  # !!!!!! chyba
    }
  }
  return(bestMotif)
}

