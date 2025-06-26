# Structure tensors

*Structure tensors* ${\bf a}^{(k)}$ are defined as the mean of the $k$-th repeated outer product of a crystallographic axis with itself, say $\bf n$:

$$ 
{\bf a}^{(k)} := \langle {\bf n}^k \rangle.
$$

In the case of a discrete ensemble of $N$ grains with slip plane normals ${\bf n}_i$ for $i=1,\cdots,N$, the structure tensors are 

$$ 
{\bf a}^{(k)} = \dfrac{\sum_i w_i {\bf n}_i^k}{\sum_i w_i},
$$

where $w_i$ is the *weight* attributed to the $i$-th grain, typically its size. 

Alternatively, if the distribution function $n(\theta,\phi)$ of ${\bf n}$-axes is known, the structure tensors are 

$$ 
{\bf a}^{(k)} = \dfrac{\int_{S^2} n \hat{{\bf r}}^k \, \mathrm{d}\Omega}{\int_{S^2} n \, \mathrm{d}\Omega} 
,
$$

where $\mathrm{d}\Omega = \sin(\theta) \mathrm{d}\theta \mathrm{d}\phi$ is the infinitesimal solid angle and $\hat{{\bf r}}(\theta,\phi)$ is the radial unit vector.

## Principal frame

Since $n(\theta,\phi)$ and $b(\theta,\phi)$ are antipodally symmetric, odd moments (odd $k$) vanish identically. Hence, ${\bf a}^{(2)}(n)$ and  ${\bf a}^{(2)}(b)$ measure the variance of $\bf n$- and $\bf b$-axes, respectively, around the three coordinate axes.

Posing ${\bf a}^{(2)}$ in its principal frame

$$ 
{\bf a}^{(2)} = \left[\begin{matrix}
\lambda_1 & 0 & 0\\ 
0 & \lambda_2 & 0\\ 
0 & 0 & \lambda_3\\ 
\end{matrix}\right]
$$
    
therefore has a similar interpretation as in PCA:
the first principal component (eigenvector ${\bf m}_1$) is the direction that maximizes the variance (eigenvalue $\lambda_1$) of the projected data (red curve), the second component is the direction orthogonal to the first component that maximizes the variance of the projected data, and so on with the third component. 

/// html | div[style='width: 100%; text-align: center;']
![](https://raw.githubusercontent.com/nicholasmr/specfab/refs/heads/main/images/harmonic-expansion/a2.png#center){: style="width:600px"}
///

## Convert to spectral 

Converting between spectral and tensorial representations is a linear problem in the sense that 

$$
{\bf a}^{(k)} = {\bf f}(\hat{n}_2^{m}, \hat{n}_4^{m}, \cdots, \hat{n}_k^{m}) 
,
\qquad\text{(for all $m$)}
$$

where ${\bf f}$ is linear in its arguments and $\hat{n}_l^m = n_l^m/n_0^0$.

In the case of ${\bf a}^{(2)}$, the relation is simple:

$$
{\bf a}^{(2)} = \frac{{\bf I}}{3} + \sqrt{\frac{2}{15}}
\left[\begin{matrix}
\operatorname{Re}[\hat{n}_2^2] - \dfrac{1}{2}\sqrt{\dfrac{2}{3}} \hat{n}_2^0 & -\operatorname{Im}[\hat{n}_2^2] & -\operatorname{Re}[\hat{n}_2^1] \\ 
 & -\operatorname{Re}[\hat{n}_2^2] - \dfrac{1}{2}\sqrt{\dfrac{2}{3}} \hat{n}_2^0  & \operatorname{Im}[\hat{n}_2^1] \\ 
\mathrm{sym.} &  & \sqrt{\dfrac{2}{3}} \hat{n}_2^0
\end{matrix}\right]
,
$$

but for higher-order structure tensors the expressions are long (not shown).
The above applies to $b(\theta,\phi)$ as well.

The following code example shows how to convert between the two CPO representations:

```python
--8<-- "docs/snippets/convert-representation.py"
```

## Construct from measurements

The harmonic expansion coefficients of any CPO can be determined from discrete measurements of its crystallographic axes.
This requires constructing the corresponding structure tensors (of each crystallographic axis), from which the expansion coefficients may be derived.

In the case of glacier ice where ${\bf n} = {\bf c}$, this can be done as follows:

```python
--8<-- "docs/snippets/discrete-to-spectral.py"
```

Note that constructing ${\bf a}^{(6)}$ is to be preferred over ${\bf a}^{(4)}$ and ${\bf a}^{(2)}$, since it contains more information on the fine-scale structure of the distribution: ${\bf a}^{(6)}$ encodes information about the expansion coefficients $l\leq 6$, ${\bf a}^{(4)}$ about $l\leq 4$, and ${\bf a}^{(2)}$ about $l\leq 2$. 

