# Seedfinder

A finder for Minecraft seeds fulfilling certain criteria.

Required libraries are pthread (for multithreading) and cubiomes (for Minecraft world generation).  
Example compilation:  
`cc -Ofast -lpthread -I path/to/cubiomes/ find_base_loc_multi.c path/to/cubiomes/libcubiomes.a`
