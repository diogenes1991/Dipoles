#procedure replaceIndices(dollar, index)
#ifdef `dollar' 
    .sort
    argument;
       id `dollar' = `index'; 
    endargument;
    .sort
#endif
#endprocedure
