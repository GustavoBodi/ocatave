function result = metodos()
  n = input("");
  P = read_poly(n);
  a_b = compute_best_quota(P);
  printf("Best quotas = [%.2f, %.2f]\n", a_b(1), a_b(2));
  x = compute_root(P, a_b(1), a_b(2));
  printf("%6f\n", x);
 endfunction

# Leitura dos polinomios
function P = read_poly(n)
  P = zeros(1, n+1);
  for k = 1:(n+1)
    expPower = n - (k-1);
    P(k) = input("");
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
  r1 = quota_minima(P);
  r2 = quota_maxima(P);
  a_b = [r1, min([upper_kojima, r2, upper_cauchy])];
endfunction

# Método de Cauchy
function result = cauchy(P)
  item_max = 5000;
  tol = 1e-6;

  n = length(P) - 1;
  a1 = P(1);

  x = 0;
  for k = 1:item_max
    soma = 0;
    for i = 2:length(P)
      expoente = n - i + 1;
      soma += abs(P(i) / a1) * power(x, expoente);
    endfor
    x1 = power(soma, (1 / n));
    if (abs(x1 - x) < tol)
      break;
    endif
    x = x1;
  endfor

  result = x;
endfunction

# Método de Kojima
function result = kojima(P)
  a1 = P(1);
  n = length(P);
  resultados = zeros(1,n);
  for i = 1:(n-1)
    ai1 = P(i+1);
    resultados(i) = power(abs(ai1 / a1), (1/i));
  endfor
  sorted = fliplr(sort(resultados));
  result = sorted(1) + sorted(2);
 endfunction

 # Cálculo da quota minima
function result = quota_minima(P)
  den = (max(abs(P(1:end-1))) / abs(P(end))) + 1;
  result = 1 / den;
endfunction

 # Cálculo da quota maxima
function result = quota_maxima(P)
  max_coef = max(abs(P(2:end))) / abs(P(1)) + 1;
  result = max_coef;
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

