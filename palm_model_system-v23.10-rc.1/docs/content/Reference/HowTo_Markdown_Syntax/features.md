# Documentation Features

this is a subtitle.

---

## Text and Images

![Placeholder](https://dummyimage.com/300x200/eee/aaa){ align=right }

Lorem ipsum dolor sit amet, consectetur adipiscing elit. ^^Nulla et euismod nulla.^^ Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. ^^Nulla et euismod nulla.^^ Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa.

## Links

- Open our old documentation [here](https://palm.muk.uni-hannover.de/trac) in the same tab.
- Open our old documentation [here](https://palm.muk.uni-hannover.de/trac){ target=_blank } in a new tab.

## Buttons

[trac wiki](https://palm.muk.uni-hannover.de/trac){ .btn .btn-primary target=_blank }

## LaTeX Equations

$$
p(x|y) = \frac{p(y|x)p(x)}{p(y)}, \left(p(x|y) = \frac{p(y|x)p(x)}{p(y)}\right).
$$

$$
E(\mathbf{v}, \mathbf{h}) = -\sum_{i,j}w_{ij}v_i h_j - \sum_i b_i v_i - \sum_j c_j h_j
$$

this is text with inline LaTeX math $\left[3 < 4 \right]$ and $\frac{p(y|x)p(x)}{p(y)}$ and so on...

\begin{align}
    p(v_i=1|\mathbf{h}) & = \sigma\left(\sum_j w_{ij}h_j + b_i\right) + \frac{p(y|x)p(x)}{p(y)} \\
    p(h_j=1|\mathbf{v}) & = \sigma\left(\sum_i w_{ij}v_i + c_j\right)
\end{align}

## Code

``` fortran
FUNCTION magnus_0d( t )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) ::  t  !< temperature (K)

    REAL(wp) ::  magnus_0d

!
!-- Saturation vapor pressure for a specific temperature:
    magnus_0d =  611.2_wp * EXP( 17.62_wp * ( t - degc_to_k ) / ( t - 29.65_wp  ) )

 END FUNCTION magnus_0d
```

``` python
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

The `#!python range()` function is used to generate a sequence of numbers.


``` c
#include <stdio.h>

int main(void) {
  printf("Hello world!\n");
  return 0;
}
```

``` c++
#include <iostream>

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit[^1].

int main(void) {
  std::cout << "Hello world!" << std::endl;
  return 0;
}
```

## Text Formating

![Placeholder](https://dummyimage.com/300x200/eee/aaa){ align=left }

Text can be {--deleted--} and replacement text {++added++}. This can also be
combined into {~~one~>a single~~} operation. {==Highlighting==} is also
possible {>>and comments can be added inline<<}.

Superscript A^T^A and subscript H~2~0 are also possible.

Lorem ipsum dolor sit amet, consectetur adipiscing elit. ^^Nulla et euismod nulla.^^ Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa.

## Tables

| Method      | Description      |
| :---------- | :--------------- |
| `GET`       | Fetch resource   |
| `PUT`       | Update resource  |
| `DELETE`    | Delete resource  |


## Abbreviations

The HTML specification is maintained by the W3C.

*[HTML]: Hyper Text Markup Language
*[W3C]: World Wide Web Consortium

## Lists

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit[^1].



* [x] Lorem ipsum dolor sit amet, consectetur adipiscing elit
* [x] Nulla lobortis egestas semper
* [x] Curabitur elit nibh, euismod et ullamcorper at, iaculis feugiat est
* [ ] Vestibulum convallis sit amet nisi a tincidunt
    * [x] In hac habitasse platea dictumst
    * [x] In scelerisque nibh non dolor mollis congue sed et metus
    * [x] Sed egestas felis quis elit dapibus, ac aliquet turpis mattis
    * [ ] Praesent sed risus massa
* [ ] Aenean pretium efficitur erat, donec pharetra, ligula non scelerisque
* [ ] Nulla vel eros venenatis, imperdiet enim id, faucibus nisi



## Admonitions

!!! note
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa.

!!! warning
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa.

!!! danger
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor massa, nec semper lorem quam in massa.


## How to include other files

{{ include_markdown('missing_file.md') }}


[^1]:
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla et euismod
    nulla. Curabitur feugiat, tortor non consequat finibus, justo purus auctor
    massa, nec semper lorem quam in massa.
