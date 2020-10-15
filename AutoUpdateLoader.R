#This funciton will call library, and update all required packages
#AULibrary == Auto Updater Library call
AULibrary <- function(expr){
  
  # Function to capture error text - included inline to keep it out of other namespaces
  # https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
  catchToList <- function(expr) {
    val <- NULL
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, w$message)
      invokeRestart("muffleWarning")
    }
    myError <- NULL
    eHandler <- function(e) {
      myError <<- e$message
      NULL
    }
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
    list(value = val, warnings = myWarnings, error=myError)
  } 
  ###################################################################################
  
  # Get the error
  error <- catchToList(expr)$error
  
  #If it does error parse for the "installed before R 4.0.0" problem
  if(!is.null(error)){
    
    # If it throws a different error - print the error and quit
    if(!grepl('installed before R 4.0.0', error)){
      print(paste('Unsupported Error:', error))
      return()
    }
    
    # While it keep throwing the error install each package and retry loading
    while(grepl('installed before R 4.0.0', error)){
      pkg <- gsub(".* package (.*) was installed.*", '\\1', error)
      pkg <- gsub("(^.)|(.$)", '', pkg)
      install.packages(pkg)
      
      error <- catchToList(expr)$error
      
      # If there are no more errors, quit
      if(is.null(error)){
        return()
      }
    }
  }
}
