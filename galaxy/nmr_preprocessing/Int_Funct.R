beginTreatment <- function(name, Signal_data = NULL, Signal_info = NULL, 
                           force.real = FALSE) {
  
  cat("Begin", name, "\n")
  
  
  # Formatting the Signal_data and Signal_info -----------------------
  
  vec <- is.vector(Signal_data)
  if (vec) {
    Signal_data <- vec2mat(Signal_data)
  }
  if (is.vector(Signal_info)) {
    Signal_info <- vec2mat(Signal_info)
  }
  if (!is.null(Signal_data)) {
    if (!is.matrix(Signal_data)) {
      stop("Signal_data is not a matrix.")
    }
    if (!is.complex(Signal_data) && !is.numeric(Signal_data)) {
      stop("Signal_data contains non-numerical values.")
    }
  }
  if (!is.null(Signal_info) && !is.matrix(Signal_info)) {
    stop("Signal_info is not a matrix.")
  }
  
  
  Original_data <- Signal_data
  
  # Extract the real part of the spectrum ---------------------------
  
  if (force.real) {
    if (is.complex(Signal_data)) {
      Signal_data <- Re(Signal_data)
    } else {
      # The signal is numeric Im(Signal_data) is zero anyway so let's avoid
      # using complex(real=...,imaginary=0) which would give a complex signal
      # in endTreatment()
      force.real <- FALSE
    }
  }
  
  
  # Return the formatted data and metadata entries --------------------
  
  return(list(start = proc.time(), vec = vec, force.real = force.real, 
              Original_data = Original_data, Signal_data = Signal_data, Signal_info = Signal_info))
}

binarySearch <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target 
  # if it doesn't exists and lower, it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target 
  # a should be monotone but can be increasing or decreasing
  
  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  N <- length(a)
  if ((a[N] - target) * (a[N] - a[1]) <= 0) {
    return(N)
  }
  if ((a[1] - target) * (a[N] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- N
  while (amin + 1 < amax) {
    amid <- floor((amin + amax)/2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # a[amid] == a[amax] or a[amid] == target In both cases, a[amid] ==
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # (lower && a[amin] < a[amax]) || (!lower && a[min] > a[max]) 
    # If increasing and we want the lower, we take amin 
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}

checkArg <- function(arg, checks, can.be.null=FALSE) {
  check.list <- list(bool=c(is.logical, "a boolean"),
                     int =c(function(x){x%%1==0}, "an integer"),
                     num =c(is.numeric, "a numeric"),
                     str =c(is.character, "a string"),
                     pos =c(function(x){x>0}, "positive"),
                     pos0=c(function(x){x>=0}, "positive or zero"),
                     l1 =c(function(x){length(x)==1}, "of length 1")
  )
  if (is.null(arg)) {
    if (!can.be.null) {
      stop(deparse(substitute(arg)), " is null.")
    }
  } else {
    if (is.matrix(arg)) {
      stop(deparse(substitute(arg)), " is not scalar.")
    }
    for (c in checks) {
      if (!check.list[[c]][[1]](arg)) {
        stop(deparse(substitute(arg)), " is not ", check.list[[c]][[2]], ".")
      }
    }
  }
}

endTreatment <- function(name, begin_info, Signal_data) {
  
  # begin_info: object outputted from beginTreatment
  
  
  # Formatting the entries and printing process time -----------------------
  end_time <- proc.time()  # record it as soon as possible
  start_time <- begin_info[["start"]]
  delta_time <- end_time - start_time
  delta <- delta_time[]
  cat("End", name, "\n")
  cat("It lasted", round(delta["user.self"], 3), "s user time,", round(delta["sys.self"],3),
      "s system time and", round(delta["elapsed"], 3), "s elapsed time.\n")
  
  
  if (begin_info[["force.real"]]) {
    # The imaginary part is left untouched
    i <- complex(real = 0, imaginary = 1)
    Signal_data <- Signal_data + i * Im(begin_info[["Original_data"]])
  }
  
  if (begin_info[["vec"]]) {
    Signal_data <- Signal_data[1, ]
  }
  
  # Return the formatted data and metadata entries --------------------
  return(Signal_data)
}

Interpol <- function(t, y) {
  # y: sample
  # t : warping function
  
  m <- length(y)
  # t <= m-1
  # because if t > m-1, y[ti+1] will be NA when we compute g
  valid <- 1 <= t & t <= m-1 # FIXME it was '<' in Bubble v2
  s <- (1:m)[valid]
  ti <- floor(t[s])
  tr <- t[s] - ti
  g <- y[ti + 1] - y[ti]
  f <- y[ti] + tr * g
  list(f=f, s=s, g=g)
}

vec2mat <- function(vec) {
  
  return(matrix(vec, nrow = 1, dimnames = list(c(1), names(vec)))) 
  
}

getArg <- function(arg, info, argname, can.be.absent=FALSE) {
  if (is.null(arg)) {
    start <- paste("impossible to get argument", argname, "it was not given directly and");
    if (!is.matrix(info)) {
      stop(paste(start, "the info matrix was not given"))
    }
    if (!(argname %in% colnames(info))) {
      if (can.be.absent) {
        return(NULL)
      } else {
        stop(paste(start, "is not in the info matrix"))
      }
    }
    if (nrow(info) < 1) {
      stop(paste(start, "the info matrix has no row"))
    }
    arg <- info[1,argname]
    if (is.na(arg)) {
      stop(paste(start, "it is NA in the info matrix"))
    }
  }
  return(arg)
}

binarySearch <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target 
  # if it doesn't exists and lower, it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target 
  # a should be monotone but can be increasing or decreasing
  
  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  N <- length(a)
  if ((a[N] - target) * (a[N] - a[1]) <= 0) {
    return(N)
  }
  if ((a[1] - target) * (a[N] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- N
  while (amin + 1 < amax) {
    amid <- floor((amin + amax)/2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # a[amid] == a[amax] or a[amid] == target In both cases, a[amid] ==
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # (lower && a[amin] < a[amax]) || (!lower && a[min] > a[max]) 
    # If increasing and we want the lower, we take amin 
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}

# returns the discrete borders of the interval for a numeric vector a

indexInterval <- function (a, from, to, inclusive=TRUE) {
  # If inclusive and from <= to, we need to take the lower
  # If not inclusive and from > to, we need to take the lower too
  lowerFrom <- xor(inclusive, from > to)
  fromIndex <- binarySearch(a, from, lowerFrom)
  toIndex <- binarySearch(a, to, !lowerFrom)
  return(fromIndex:toIndex)
}
