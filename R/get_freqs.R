.get_freqs <- function(data) {
    if (sum(as.numeric(data))==1) {
        warning('I assume you provided me a list of relative frequences')
        return(data)
    } else {
        warning(paste('I assume you provided me a list of absolute values.\n',
                      'I internally transformed them in relative frequencies.\n',
                      'Otherwise check in the original data',
                      ' why their sum is not 1.'))
        freqs <- data/sum(as.numeric(data))
    }
    return(freqs)
}
