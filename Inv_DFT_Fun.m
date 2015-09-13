function Inv_DFT_Image=Inv_DFT_Fun(Image)

i=sqrt(-1);
Inv_DFT_Image=Image;

[Row_num , Column_num]=size(Image);
for ii=1:Row_num
    for ij=1:Column_num
        sum_n=0;
        for fur_i=1:Column_num
            sum_n=sum_n+Image(ii,fur_i)*exp(2*pi*i*(ij-1)*(fur_i-1)/Column_num);
        end
        Inv_DFT_Image(ii,ij)=1/sqrt(Column_num)*sum_n;
    end
end