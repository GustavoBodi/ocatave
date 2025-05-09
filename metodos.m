function result = metodos()
  n = input("");
  disp(n);
  P = read_poly(n);
  a_b = compute_best_quota(P);
  x = compute_root(P, a_b(1), a_b(2));
  printf("%6f\n", x);
 endfunction

# Leitura dos polinomios
function P = read_poly(n)
  P = zeros(1, n+1);
  for k = 1:(n+1)
    expPower = n - (k-1);
    prompt = sprintf('Insira o coeficiente para a%d (for x^%d): ', k, expPower);
    P(k) = input(prompt);
  end
endfunction

# Escolha dos procedimentos para computacao final
function root = compute_root(P, a, b)
  if (computar_func(P, a) * computar_func(P, b)) <= 0
    root = metodo_falsa_posicao(P, a, b);
  else
    root = metodo_secantes(P, a, b);
  endif
endfunction

# Computar o minimo e o máximo intervalos aproximados
function a_b = compute_best_quota(P);
  upper_kojima = kojima(P);
  upper_cauchy = cauchy(P);
  candidata_minima = quota_minima(P);
  candidata_maxima = quota_maxima(P);
  a_b = [max(candidata_maxima, candidata_minima), min(upper_kojima, upper_cauchy)];
endfunction

# Método de Cauchy
function result = cauchy(P)
  an = P(1);
  lista_divisao = P(2:end) / an;
  result = 1 + max(lista_divisao);
endfunction

# Método de Kojima
function result = kojima(P)
  coeficientes = P;
  an = coeficientes(1);
  n = length(coeficientes);
  resultados = zeros(1,n);
  for k = 1:(n)
    ak = coeficientes(k);
    resultados(k) = power(abs(ak / an), (1/k));
  endfor
  result = max(resultados);
 endfunction

 # Cálculo da quota minima
 function result = quota_minima(P)
   coeficientes = P;
   a0 = coeficientes(end);
   n = length(P) - 1;
   raizes_candidatas = zeros(1,n);
   for k = 1:(n)
     ak = coeficientes(k);
     raizes_candidatas(k) = power(abs( ak / a0 ), (1 / (n - k)));
   endfor
   minimo = min(raizes_candidatas);
   result = minimo / (1 + minimo);
endfunction

 # Cálculo da quota maxima
function result = quota_maxima(P)
  coeficiente_invertidos = fliplr(P);
  r_rev = cauchy(coeficiente_invertidos);
  result = 1 / r_rev;
endfunction

# Aplicacao secantes
function resultado = metodo_secantes(P, a, b)
  x0 = (a + b) / 2;
  x1 = x0 + 0.01;
  x2 = 0;
  iter_max = 5000;
  for j = 1:(iter_max)
    x2 = x1 - computar_func(P, x1) * (x1 - x0) / (computar_func(P, x1) - computar_func(P, x0));
    x0 = x1;
    x1 = x2;
    if (abs(x0 - x1) < power(10, -6))
      break
    endif
  endfor
  resultado = x2;
endfunction

# Aplicacao falsa posicao
function raiz = metodo_falsa_posicao(P, a, b)
  tol = 1e-6;
  iter_max = 5000;
  fa = computar_func(P,a);
  fb = computar_func(P,b);
  for i = 1:iter_max
    c = (a*fb - b*fa)/(fb - fa);
    fc = computar_func(P,c);
    if abs(fc) < tol
      break;
    endif
    if fa*fc < 0
      b = c; fb = fc;
    else
      a = c; fa = fc;
    endif
  endfor
  raiz = c;
endfunction

# Computar o valor do polinomio para x
function y = computar_func(P, x)
  y = P(1);
  n = length(P);
  for k = 2:n
    y = y * x + P(k);
  endfor
endfunction

