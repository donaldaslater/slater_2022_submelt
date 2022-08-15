EN4 is an objective ocean analysis based on profile data. I used version 4.2.1 of the objective analysis with Gouretski and Reseghetti (2010) corrections. Full description at https://www.metoffice.gov.uk/hadobs/en4/index.html

To download the data run this on the command line  
`wget -i EN4_downloadlist.txt`  
It's worth checking this has worked as sometimes it can download only half the file and it's hard to notice. Unpack the files using  
`unzip "*.zip"`  
Once unpacked, the EN4 monthly files for 1979-2020 occupy ~30GB.

Running the `get_oceanpoints.m` script just requires one of the monthly EN4 files, e.g. EN.4.2.1.f.analysis.g10.202006.nc
