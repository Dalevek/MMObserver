%% Generowanie wartoœcio pocz¹tkowych x0 %%
x0 = eye(size(A,1));
x0_d = zeros(size(A,1),size(A,1)+1);
x0_d(1:size(A,1),1:size(A,1)) = x0;
%% Wygenerowanie wektora x0 do symulacji
x0=reshape(x0_d,[],1)

%% Generowanie macierzy A_d %%
s = (size(A,1)+1)*size(A);
A_d = zeros(s);
n = size(A,1);
for count = 1:n+1
    n_s = (count-1)*n+1;
    n_e = n_s+n-1;
    A_d(n_s:n_e,n_s:n_e) = A;
end
%% Zapis wygenerowanej macierzy do macierzy A %%
A=A_d

%% Generowanie macierzy B_d %%
s = (size(B,1)+1)*size(B); 
B_d = zeros(s);

row = size(B_d,1);
col = size(B_d,2);
row_B = size(B,1);
col_B = size(B,2);

col_temp=1;
for r = 1:row_B:row
    B_d(r:(row_B+r-1),col_temp:(col_B+col_temp-1)) = B;
    col_temp=col_temp+col_B;
end
%% Zapis wygenerowanej macierzy do macierzy B %%
B=B_d



%% Generowanie macierzy C_d %%
s = (size(C,1))*(1+size(C,2));
C_d = zeros(s);

row = size(C_d,1);
col = size(C_d,2);
row_C = size(C,1);
col_C = size(C,2);

col_temp=1;
for r = 1:row_C:row
    C_d(r:(row_C+r-1),col_temp:(col_C+col_temp-1)) = C;
    col_temp=col_temp+col_C;
end
%% Zapis wygenerowanej macierzy do macierzy C %%
C=C_d


%% Generowanie macierzy L_d %%
s = (size(L,1)+1)*size(L);
L_d = zeros(s);

row = size(L_d,1);
col = size(L_d,2);

row_L = size(L,1);
col_L = size(L,2);

col_temp=1;
for r = 1:row_L:row
    L_d(r:(row_L+r-1),col_temp:col_L+col_temp-1) = L;
    col_temp=col_temp+col_L;
end
%% Zapis wygenerowanej macierzy do macierzy L %%
L=L_d