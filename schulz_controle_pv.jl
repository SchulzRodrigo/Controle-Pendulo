##################################
#   Autor: Rodrigo André Schulz  #
##################################

using Plots
using LinearAlgebra
using XLSX
##################################################################
function malha(ti,tf)
    #malha
    global h=0.001 #tamanho do passo
    n=round(Int32,(tf-ti)/h) #divide o intervalo de tempo total pelo tamanho do passo h
    n=n+1
    return n,h
end
##################################################################
function RK4(ti,tf,y0,F) #método Runge-Kutta de ordem 4
    (n,h)=malha(ti,tf)
    t=zeros(n)
    y=zeros(n,length(y0))
    y[1,:]=y0
    t=[ti:h:tf;]
    for i=2:n
        K1=h*F(t[i-1],y[i-1,:])
        K2=h*F(t[i-1]+h/2,y[i-1,:]+K1/2)
        K3=h*F(t[i-1]+h/2,y[i-1,:]+K2/2)
        K4=h*F(t[i-1]+h,y[i-1,:]+K3)
        y[i,:]=y[i-1,:]+(1.0/6.0)*(K1+(2*K2)+(2*K3)+K4)
    end
    return (t,y)
end
#########################################################
####################################################################
function int_simpson_comp(f,n,a,b) #integral de simpson
    if n%2==1 #se n é ímpar
        n=n+1
        println("n dever ser par. Utilizamos n=$n")
    end
    h=(b-a)/n
    int=f(a)+f(b)
    x=a+h
    for i=1:n-1
        if i%2==1 #se i é ímpar
            int+=4*f(x)
        else #se i é par
            int+=2*f(x)
        end
        x+=h
    end
    return int=int*h/3
end
#######################################################
function L(t) #Cúbica ou reta
    if t<=ti
        aux=Li
    elseif ti<t<tf
        aux = coef_cubica[1]*t^3 + coef_cubica[2]*t^2 + coef_cubica[3]*t + coef_cubica[4] # Cúbica
        #aux=((Lf-Li)/(tf-ti))*(t-ti)+Li #equação da reta
    else
        aux = Lf
    end
    return aux
end
#######################################################
#X'(t)=AX(t)+Bu(t); X(ti)=X0;
# Pendulo + 
ti=0
tf=4
Xi=[0;0;0;0]
Xf=[2;0;2;0]  
Li=1.0
Lf=0.5
#Coeficientes Cúbica
Matriz_Coeficientes = [3*ti^2 2*ti 1 0; 3*tf^2 2*tf 1 0; ti^3 ti^2 ti 1; tf^3 tf^2 tf 1]
Matriz_independentes = [0;0;Li;Lf]
coef_cubica = Matriz_Coeficientes\Matriz_independentes
A=[0.0 1 0 0;0 0 0 0;0 0 0 1;10/L(ti) 0 -10/L(ti) 0] 
B=[0; 0.1; 0; 0] 
#######################################################
#Grafico L(t)
pyplot()
theme(:default)
grafico_L = plot(L, 
    title = "Variação do Comprimento da Haste com o Tempo",
    xlabel = "tempo t",
    ylabel = "L(t)",
    leg = false,
    xlim = (ti, tf)
)
savefig(grafico_L, "grafico_L_pv.png")
#######################################################
function At(t)
    aux=A
    aux[4,:]=[10/L(t) 0 -10/L(t) 0]
    return aux
end
############################################
function u(s)
    if 0<s<=tf
        aux=adjoint(B)*exp(adjoint(ArgA)*(tf-s))*inv(ArgW_T)*(Xf-exp(ArgA*tf)*Xi)
    else
        aux=0
    end
    return aux
end
########################################################
function EDO(s,ArgX) #X'(t)=A(t)X(t)+Bu(t)
    X=At(s)*ArgX+B*u(s)
    return X
end
########################################################
num_testes = 100
matriz_int = zeros(num_testes, 4) #[tempo, integral posição, integral velocidade, integral soma]
for i = 1:num_testes
    println("i = ", i)
    matriz_int[i,1]=ti+i*(tf-ti)/num_testes
    global ArgA=At(matriz_int[i,1])
    M(s)=exp(ArgA*s)*B*adjoint(B)*exp(adjoint(ArgA)*s)
    global ArgW_T=int_simpson_comp(M,(tf-ti)*5,ti,tf)
    (t,X)=RK4(ti,tf+0.5*tf,Xi,EDO)
    C = findall(x->x>=tf, t) #encontra indice a partir de tf
    for j = C[1]:length(X[:,1])
        matriz_int[i,2] += abs(Xf[3]-X[j,3])*h #posição
        matriz_int[i,3] += abs(Xf[4]-X[j,4])*h #velocidade
    end
    matriz_int[i,4] = matriz_int[i,2] + matriz_int[i,3] #soma
end
matriz_int
#encontra a posição na matriz do menor valor em int_v e int_p
(min_int_v,key_min_v)= findmin(matriz_int[:, 2])
(min_int_p,key_min_p)=findmin(matriz_int[:, 3])
(min_int_s,key_min_s)=findmin(matriz_int[:, 4])
matriz_int[key_min_v,:]
matriz_int[key_min_p,:]
matriz_int[key_min_s,:]

# Substitui o tempo para fazer os gráficos de u(t) e da solução da EDO
ArgA=At(matriz_int[key_min_s,1]) # argumento é o tempo que a ser usado para calcular o comprimento da haste

M(s)=exp(ArgA*s)*B*adjoint(B)*exp(adjoint(ArgA)*s)
ArgW_T=int_simpson_comp(M,(tf-ti)*5,ti,tf)
(t,X)=RK4(ti,tf+0.5*tf,Xi,EDO)
Xt=[t X]

# Gráfico u(t)
pyplot()
theme(:default)
grafico_u = plot(u, 
    title = "Melhor Função Controle Encontrada",
    xlabel = "tempo t",
    ylabel = "u(t)",
    leg = false,
    xlim = (ti, tf+(tf-ti)/2)
)
savefig(grafico_u, "grafico_u_pv.png")

# Gráfico Pêndulo
pyplot()
theme(:default)
grafico_com_u_pv=plot(Xt[:,1], [Xt[:,2] Xt[:,3] Xt[:,4] Xt[:,5]],
title  = "Solução com menor erro",
xlabel = "tempo",
ylabel = "valor",
label  = ["posição carro" "velocidade carro" "posição objeto" "velocidade objeto"]
)
savefig(grafico_com_u_pv, "grafico_com_u_pv.png")

#Grafico Erro x Comprimento Haste
# aqui calculamos controles para diferentes valores de comprimento de haste fixada
pyplot()
theme(:default)
grafico_erro_pv=plot(map(L,matriz_int[:,1]), [matriz_int[:,2] matriz_int[:,3] matriz_int[:,4]],
title  = "Erro x Comprimento da Haste",
xlabel = "Comprimento Haste",
ylabel = "Erro Integral",
label  = ["Integral posição" "Integral velocidade" "Integral soma"]
)
savefig(grafico_erro_pv, "grafico_erro_pv.png")


