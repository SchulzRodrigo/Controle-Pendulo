##################################
#                                #
#   Autor: Rodrigo André Schulz  #
#                                #
##################################

using Plots
using LinearAlgebra #para calcular exponenciais de matrizes
using QuadGK #para calcular integrais
using DifferentialEquations #para resolver a EDO

# Dados de entrada para o sistema na forma matricial
# Pendulo + Carrinho 
#X'(t)=AX(t)+Bu(t); X(ti)=X0;
A=[0.0 1 0 0;0 0 0 0;0 0 0 1;10 0 -10 0]
B=[0; 0.1; 0; 0]
ti=0            #instante inicial
Xi=[0;0;0;0]    #dado inicial
tf=2            #tempo para controlar
Xf=[1;0;1;0]    #estado final desejado

#escrevendo a matriz do integrando de W_T
M(s)=exp(A*s)*B*adjoint(B)*exp(adjoint(A)*s)

#calculando a integral de M para obter W_T
(W_T, erro_int)=quadgk(M, ti, tf, rtol=1e-5) 


# Definindo a função de controle
function u(s)
    if 0<s<=tf
        aux=adjoint(B)*exp(adjoint(A)*(tf-s))*inv(W_T)*(Xf-exp(A*tf)*Xi)
    else
        aux=0
    end
    return aux
end

# Gerando o gráfico da função controle u(t)
grafico_u = plot(u, 
    title = "Função Controle",
    xlabel = "tempo t",
    ylabel = "u(t)",
    leg = false,
    xlim = (ti, tf+(tf-ti)/2)
)
savefig(grafico_u, "grafico_u.png")

# Definindo a EDO
function Pendulo(dX,X,p,t)
    dX[1] = X[2]
    dX[2] = p[1]*u(t)
    dX[3] = X[4]
    dX[4] = p[2]*X[3]+p[3]*X[1]
end
p=[B[2]; A[4,3];A[4,1]]
sist = ODEProblem(Pendulo,Xi,(ti,tf+(tf-ti)/2),p)

# Resolvendo a EDO e gerando o gráfico
sol = solve(sist,Tsit5())
grafico_com_u=plot(sol,
    linewidth=1,
    title = "Pêndulo Controlado",
    label=["Posição carro" "Velocidade carro" "Posição objeto" "Velocidade objeto"],
    xlabel = "tempo t",
    ylabel = "valores"
)
savefig(grafico_com_u, "grafico_com_u.png")

#Resolvendo a EDO sem o Controle
u(s)=0
sol = solve(sist,Tsit5())
grafico_sem_u=plot(sol,
    linewidth=1,
    title = "Pêndulo sem Controle",
    label=["Posição carro" "Velocidade carro" "Posição objeto" "Velocidade objeto"],
    xlabel = "tempo t",
    ylabel = "valores"
)
savefig(grafico_sem_u, "grafico_sem_u.png")
