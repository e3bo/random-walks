FROM rhub/r-minimal

RUN installr -d \
   -t "R-dev file automake autoconf linux-headers libxml2-dev" \
   -a "libxml2 icu-libs" `# needed for xml2 and stringr`  \
   tidyverse/readxl `# CRAN version (<=1.3.1) does not compile on alpine` \
   tidyverse

# above copied from github.com/r-hub/r-minimal/examples/tidyverse/Dockerfile by edwindj
#RUN apk update && apk add cmake bash # needed to build arrow library
RUN installr -d -t "cmake bash" arrow

COPY hopkins-ts-reformat.R /opt/ml/hopkins-ts-reformat.R

ENTRYPOINT ["/usr/local/bin/Rscript", "/opt/ml/hopkins-ts-reformat.R", "--no-save"]
