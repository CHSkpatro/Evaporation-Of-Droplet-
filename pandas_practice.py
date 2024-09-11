import pandas as pd
import numpy as np
data=np.arange(0,200).reshape(20,10)
print(data)
print(type(data))
#for row indices without typing all 20 manually
row=[f"Row{i}" for i in range(1,21)]#fstrings 
column=[f"Column{j}" for j in range(1,11)]
df=pd.DataFrame(data,index=row,columns=column)#making a dataframe from array 
print(type(df))
print(df)
print(df.head())#gives only top5 rows
print(df.tail())#gives only bottom5 rows
print(df.info())
print(df.describe())#gives mean std.deviation datacount percentile min max, for only those columns will appear in which int float type data is there 

#indexing using column name
print(df["Column1"])
print(type(df["Column1"]))#series, when either only one row  or column is being accessed
print(df[["Column1","Column3","Column5","Column6"]])
print(type(df[["Column1","Column3","Column5","Column6"]]))#dataframe

#using row index [Loc]
print(df.loc['Row4'])
print(df.loc[['Row2','Row4','Row8','Row10']])#these are the row names not the row indeces

#using row index column index[iloc]
print(df.head())
print(df.iloc[1:6,2:8])#will give the dataframe of 2nd row to 6th row and from 3rd column to 8th  column , indeces start from 0 and ,l:r+1
print(df.iloc[16:,:4])
print(df.iloc[:,[0,8]])#will give the entirety of column 1 and column 9 
print(df.iloc[[4,5,18,8],[3,8,5,9]])#means will give row 5 values for column4, will give row6 values for column9 and so on
#convert dataframe into arrays
print(df.iloc[:,:].values)
print(type(df.iloc[:,:].values))
print(df.iloc[[1,3,6,9],[1,2,3,4,5,6,7]].values)
print(type(df.iloc[[1,3,6,9],[1,2,3,4,5,6,7]].values))
#basic functions with df
print(df.isnull().sum())
#lets create data frame from strings
row1=[f"Row{k}" for k in range(1,3)]
col2=[f"Column{l}" for l in range(1,4)]
df2=pd.DataFrame(data=[[1,np.nan,2],[1,3,4]],index=row1,columns=col2)
print(df2)
print(df2.isnull().sum())#will give the number of null values in a column 
print(df2.info())
print(df2.isnull())#will show true if the value is null somewhere and false where the value is not null at the position of that element
print(df2.isnull().sum()==0)#will tell in which column the null value is present
#now lets go back to df dataframe
print(df["Column8"].value_counts())#will tell how many times A VALUE IS PREsent in column8
print(df["Column9"].unique())#will show the aray of unique values present in column9
print(df>2)#will show if the elements in ieach cell is greater than 2 or not 
print(df["Column2"]>2)#will show if the elements in column 2 are greater than 2 or not by true false 
print(df[df["Column5"]>49])#will print the entire row in which the 5th column element is greater than 49
print(df[(df["Column6"]>35) & (df["Column6"]<89)])#will print the entire row in dataframe where the element in the 6th column is greater than 35 and less than 89
#opening a comma separated values(csv) file witht he help of pandas
from io import StringIO
df3=pd.read_csv('industry.csv')
print(df3)
print(df3.head())
data2=('col1,col2,col3,col4\n'
      '1  ,  y   ,  z   ,   12\n'
       '11  ,  b   ,  c   ,   13\n'
       '111   ,  19  ,   23\n'
       '1111   ,  f   , np.Nan,  19   ')#here data is in stringtype
df4=pd.read_csv(StringIO(data2))#now the data has been converted into datframe
print(df4)
df5=df4.loc[0:2]#for .loc[] both the values are inclusive, means it will access the rows index 0 to 2
print(df5)
#to access the specific  columns 
df6=pd.read_csv(StringIO(data2),usecols=['col1','col4'])
print(df6)
#now we want df6 dataframe to be converted into csv file and saved as the given file name
df6.to_csv('test.csv')
#now when we open the csv file we will see that the rows are numbered twice in excel file so to have only the excel numbering of 1 2 3 4 and delete the index numbers
df6.to_csv('test2.csv',index=False)
print(df4.info())
#now let we want to specify the dataframe dtypes according to our need
df7=pd.read_csv(StringIO(data2),dtype='object')#we converted all the column data types into objects/string
print(df7.info())
#now we want to specify the 1st column as dtype int, 2nd column as dtype float, 3rd column as dtype object,and want to leave the 4 th column as it is , now when we processed it showed error as it cannot convert strings into float 
dtype={'col1':int,'col2':object,'col3':object}
df8=df4.astype(dtype)
print(df8.info())
print(df8.dtypes)#will only give us datatypes in each column
###################lets there is a link of a file , we want to convert it into csv file , we saw that the file in url is separated by tab as separator between columns 
###################pd.read_csv('link',sep='\t').to_csv('')#here csv file name we want to give it to will be here
#######################        the \t is separator symbol when tab is the separator





