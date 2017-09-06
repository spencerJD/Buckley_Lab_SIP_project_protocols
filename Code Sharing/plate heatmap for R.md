Creating a graphical heatmap of a 96-well plate:
==============================================================

>Use the following code in R to create a heatmap of a 96-well plate based on a variable like [DNA] or intensity of fluorescence

* added by Sara Sirois 9/6/2017

* taken from user "mnel" at https://stackoverflow.com/questions/13983225/plot-plate-layout-heatmap-in-r

* useful info for changing color gradient scale at http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually



This code assumes you have a dataset named "DATA" which includes the variables "ROW" (A-H), "COLUMN" (1-12), and "INTENSITY" for each sample.

### FOR CODE OUTPUT EXAMPLE, click [HERE](https://user-images.githubusercontent.com/16819535/30114617-642c1912-92e6-11e7-9437-bbeb67421194.png)

hyperlink to image: https://user-images.githubusercontent.com/16819535/30114617-642c1912-92e6-11e7-9437-bbeb67421194.png

Here's the code:

```
library(ggplot2)

ggplot(DATA, aes(y = factor(ROW, rev(levels(ROW))),x = factor(COLUMN))) + 
     geom_point(aes(colour = INTENSITY), size =9)  +  #change size of dots as appropriate
     theme_bw() + 
     labs(x=NULL, y = NULL) +   #change axis labels as appropriate
     scale_color_gradient(low="grey", high="purple")   #change color gradient scheme (unless you don't dig purple, which is incomprehensible)
     
```

