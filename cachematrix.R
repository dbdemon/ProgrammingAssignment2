#step one:create a list of 4 functions, one that sets a matrix to be a value, one that gets a matrix already set, one that sets the inverse of the matrix,
#and one that gets that inverse from the cache

makeMatrix <- function(x = matrix(numeric(0))) {
  i <- NULL
  set <- function(y) {
    x <<- y
    i <<- NULL
  }
  get <- function() x
  setinv <- function(inverse) i <<- inverse
  getinv <- function() i
  list(set = set, get = get,
       setinv = setinv,
       getinv = getinv)
}

#step two: create a function that runs the functions, first checking to see if the inverse of the matrix is already in the cache and returning that if applicable, 
#and if not, getting the matrix from the makeMatrix list, calculating the inverse, storing it in the cache, and returning it

cacheinv <- function(x) {
  i <- x$getinv()
  if(!is.null(i)) {
    message("getting cached data")
    return(i)
  }
  data <- x$get()
  i <- invert(data)
  x$setinv(i)
  i
}

#as a numeric matrix is invertable if and only if it is square and has a non-zero determinant, I decided to make a function to do the Laplace expansion of the
#matrix, remembering 2 rules: the determinant of a 1x1 matrix is the number inside it, and the determinant of an nxn matrix is the alternating dot product of 
#the entries of the first row and the determinants of their minors (the matrix you get when removing the row and column of that entry)

getdeterminant <- function(x){
  alt <- -1
  altsum <- 0
  if(length(x)==1){
    return(x)
  }
  for(y in 1:nrow(x)){
    alt <- alt*(-1)
    altsum <- altsum + alt*(x[1, y])*getdeterminant(x[-1, -1*y])
  }
  return(altsum)
}

#this function does the actual inversion of the matrix with Gaussian elimination, but first it checks the 2 qualities that a matrix needs to have to be invertible

invert <- function(x){
  #square matrices are not invertable, so if rows are not equal to columns the program helpfully reminds you of that fact
  if(nrow(x) != ncol(x)){
    message("this is not a square matrix you square, so no inverse")
    return(NULL)
  }
  #if the determinant is 0, then a matrix is not invertable, so the program checks for that first
  if(getdeterminant(x)==0){
    message("as the determinant is 0 this is a degenerate matrix you degenerate, so no inverse")
    return(NULL)
  }
  #to do a Gaussian elimination to find a matrix's inverse, I start with an identity matrix of the same size (a square matrix where if the row equals the column the
  #entry equals 1 and otherwise the entry equals 0
  invertedm <- matrix(nrow = nrow(x), ncol = ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(i == j){
        invertedm[i, j] <- 1
      }else{
        invertedm[i, j] <- 0
      }
    }
  }
  #once that's done, the Gaussian Elimination algorithm has 3 steps for every row of the matrix x and the mapped inverse inversem
  for(i in 1:nrow(x)){
    #step 1: for the ith row of matrix x, check if it has a non-zero entry in the ith column, and if not switch it and the ith row of inversem with the first row below
    #the row in x with a non-zero entry and the same row in inversem
    j <- i
    unmet <- TRUE
    while(unmet){
      if(x[i, j]==0){
        j <- j + 1
      }else{
        unmet <- FALSE
        if(i != j){
          rowx <- x[i,]
          rowinv <- invertedm[i,]
          x[i,] <- x[j,]
          invertedm[i,] <- invertedm[j,]
          x[j,] <- rowx
          invertedm[j,] <- rowinv
        }
        #step 2: divide the ith row of both matrices by what's in the ith row and column of x, making it 1
        operator <- x[i, i]
        invertedm[i,] <- invertedm[i,]/operator
        x[i,] <- x[i,]/operator
        #step 3: clear the column by subtracting the ith row for both matrices from every other row so that every other entry in the ith column of x is 0
        for(k in 1:nrow(x)){
          if(k != i){
            invertedm[k,] <- invertedm[k,] - x[k, i]*invertedm[i,]
            x[k,] <- x[k,] - x[k, i]*x[i,]
          }
        }
      }
    }
  }
  #once done, x is an identity matrix, and invertedm is the inverse of the original matrix
  return(invertedm)
}

#let's test it out: this should get the inverse of the matrix matrix(c(1, 2, 3, 4), 2, 2) by calculating it
testmatrix <- makeMatrix(matrix(c(1, 2, 3, 4), 2, 2))
cacheinv(testmatrix)
#and this should take it from the cache
cacheinv(testmatrix)
