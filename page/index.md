\enabletheorems
<!-- =============================
     ABOUT
    ============================== -->

\begin{section}{title="About this Package", name="About"}

\lead{PkgPage.jl is based upon [Franklin.jl](https://github.com/tlienart/Franklin.jl) and makes it easy to create a beautiful landing page for a package in less than 10 minutes.}

PkgPage uses [Bootstrap 4.5](https://getbootstrap.com/docs/4.5/getting-started/introduction/) and is meant to be particularly easy to use with minimal need to tweak HTML or CSS to get a great-looking result.
With it you can:

* easily insert and style figures and tables,
* easily initiate a multi-columns environment,
* easily control the layout, colors, code theme etc from the config file,
* automatic insertion of new sections into the navbar,
* automatic purging of bootstrap css to remove unused rules,
* all the goodies from [Franklin.jl](https://github.com/tlienart/Franklin.jl):
    * live-rendering,
    * support for [KaTeX](https://github.com/KaTeX/KaTeX) and [highlight.js](https://highlightjs.org/) and automatic pre-rendering,
    * and [much more](https://franklinjl.org/).

\end{section}
\begin{section}{title="The basic settings"}
Let's assume we have an underlying stock price varying as a stochastic process called $S_t$ and
let's assume without loss of generality that $S_t>0$ a.e..\\
In that case we can define as $X_t$ the following:

$$X_t=\log(S_t) $$

Let's assume we want to price and compute the sensitivities of a european option with maturity $T$ and payoff $\phi$.

We can express the price of such option as:

$$V_0=\mathbb{E}(\phi(S_T)) $$

or, in function of $X_T$:
$$V_0=\mathbb{E}(\Phi(X_T)) $$

Where:

$$\Phi(x)=\phi(e^x) $$
\end{section}
\begin{section}{title="Pathwise Method"}
Since:
$$V_0=\mathbb{E}(\Phi(X_T))$$

$$\partial_\theta V_0=\partial_\theta \mathbb{E}(\Phi(X_T)) = \partial_\theta \mathbb{E}(\dfrac{\partial \Phi(X_T)}{\partial X_T}\dfrac{\partial X_T}{\partial \theta})$$

Where the following is called tangent process:

$$Y_t=\dfrac{\partial X_t}{\partial \theta}$$

Unfortunately this method does not provide usefull results for binary options or for n-th order sensitivities, indeed in case of binary options we have:

$$\Phi(x)=\mathbb{1}(x)$$

Hence:

$$\partial_x \Phi(x)=0$$

\end{section}
\begin{section}{title="Likelihood Ratio Method"}
Let's assume we have an underlying stock price varying as a stochastic process called $S_t$ and
let's assume without loss of generality that $S_t>0$ a.e..\\
In that case we can define as $X_t$ the following:

$$X_t=\log(S_t) $$

Let's assume we want to price and compute the sensitivities of a european option with maturity $T$ and payoff $\phi$.

We can express the price of such option as:

$$V_0=\mathbb{E}(\phi(S_T)) $$

or, in function of $X_T$:
$$V_0=\mathbb{E}(\Phi(X_T)) $$

Where:

$$\Phi(x)=\phi(e^x) $$

If we develop the integral:

$$V_0=\mathbb{E}(\Phi(X_T))=\int_{\mathbb{R}}  \Phi(x) f(x) dx $$

where the following is the density of the log yield:

$$f(x)=f_{X_T}(x)$$

If we differentiate for a generic parameter $\theta$:

$$\partial_\theta V_0=\partial_\theta \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial_\theta f(x,\theta) dx$$

or:

$$\partial_\theta V_0 = \int_{\mathbb{R}}  \Phi(x) \partial_\theta(\log(f(x,\theta))) f(x,\theta) dx$$

or:

$$\partial_\theta V_0 = \mathbb{E}(\Phi(X_T) \partial_\theta \log(f(X_T,\theta)))$$

In case of n-th order sensitivities:

We observe that for a generic positive differentiable function f we can have:

$$ g(x,\theta) =  \log{f(x,\theta)}$$

or:

$$ f(x,\theta) =  e^{g(x,\theta)}$$

Then, according to [Faa di Bruno's formula](https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula):

$$\partial_{\theta^n} e^{g(x,\theta)} = e^{g(x,\theta)} B_n(\{\partial_{\theta^i} g(x,\theta)\}_{i = 1 : n})  $$

Where $B_n$ is the complete exponential [Bell polynomial](https://en.wikipedia.org/wiki/Bell_polynomials) of n-th order.
$$B_0 =   1 $$
$$B_1(x_1) =   x_1 $$
$$B_2(x_1,x_2) =   x_1^2 + x_2 $$
$$B_3(x_1,x_2,x_3) =   x_1^3 + 3x_1 x_2 + x_3 $$
$$B_4(x_1,x_2,x_3,x_4) =   x_1^4 + 6 x_1^2 x_2 + 4 x_1 x_3 + 3 x_2^2 + x_4 $$
$$B_5(x_1,x_2,x_3,x_4,x_5) =   x_1^5 + 10 x_2 x_1^3 + 15 x_2^2 x_1 + 10 x_3 x_1^2 + 10 x_3 x_2 + 5 x_4 x_1 + x_5 $$
$$B_6(x_1,x_2,x_3,x_4,x_5,x_6) =   x_1^6 + 15 x_2 x_1^4 + 20 x_3 x_1^3 + 45 x_2^2 x_1^2 + 15 x_2^3 + 60 x_3 x_2 x_1 + 15 x_4 x_1^2 + 10 x_3^2 + 15 x_4 x_2 + 6 x_5 x_1 + x_6 $$
Hence:

$$\partial_{\theta^n} f(x,\theta) = f(x,\theta) B_n(\{\partial_{\theta^i} \log{f(x,\theta)}\}_{i = 1 : n})  $$

By applying this to the original problem:

$$\partial_{\theta^n} V_0=\partial_{\theta^n} \int_{\mathbb{R}}  \Phi(x) f(x,\theta) dx = \int_{\mathbb{R}}  \Phi(x) \partial_{\theta^n} f(x,\theta) dx$$

By applying the previous result:

$$\partial_{\theta^n} V_0= \int_{\mathbb{R}}  \Phi(x) f(x,\theta) B_n(\{\partial_{\theta^i} \log{f(x,\theta)}\}_{i = 1 : n}) dx$$

or, in terms of expectation:


$$\partial_{\theta^n} V_0 = \mathbb{E}(\Phi(X_T) B_n(\{\partial_{\theta^i} \log{f(X_T,\theta)}\}_{i = 1 : n}))$$



\end{section}


\begin{section}{title="Likelihood Ratio Method for Black and Scholes"}
In case of trivial models the density is known in analytic form, otherwise one need to compute it numerically from the characteristic function.
In case of the Black and Scholes Model, and a European Option we have the following:

$$X_t=\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) t + \sigma W_t$$

$$\log(f(x,S_0,r,d,\sigma,T)) = -\dfrac{1}{2}((\dfrac{x-(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T)}{\sigma \sqrt{T}})^2 + \log(2 \pi)) - \log(\sigma\sqrt{T})$$

From this we can compute the various partial derivatives:

$$ \partial_{S_0} \log(f(x,S_0,r,d,\sigma,T)) = \frac{x - \log\left( S_0 \right) - \left(r - d - \frac{\sigma^{2}}{2} \right) T}{\sigma^{2} T S_0}$$
$$ \partial_{r} \log(f(x,S_0,r,d,\sigma,T)) = \frac{ x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) }{\sigma^{2}}$$
$$ \partial_{d} \log(f(x,S_0,r,d,\sigma,T)) = - \partial_{r} \log(f(x,S_0,r,d,\sigma,T))$$
$$ \partial_{\sigma} \log(f(x,S_0,r,d,\sigma,T)) =\frac{ - \left( x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) \right) \left( \frac{T}{\sqrt{T}} - \frac{x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right)}{\sigma^{2} T} \sqrt{T} \right)}{\sqrt{T} \sigma} + \frac{-1}{\sigma}$$
$$ \partial_{T} \log(f(x,S_0,r,d,\sigma,T)) =\frac{ - \left( x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right) \right) \left( \frac{d - r + \frac{\sigma^{2}}{2}}{\sqrt{T} \sigma} + \frac{ - \frac{x - \log\left( S_0 \right) - T \left(  r - d - \frac{\sigma^{2}}{2} \right)}{\sigma^{2} T} \sigma}{2 \sqrt{T}} \right)}{\sqrt{T} \sigma} + \frac{-1}{2 T}$$

Hence, a delta for a generic european option can be computed as the following expectation:

$$\Delta = \partial_{S_0} V_0 = \mathbb{E}(\Phi(X_T) \frac{X_T - \log\left( S_0 \right) - \left(r - d - \frac{\sigma^{2}}{2} \right) T}{\sigma^{2} T S_0} )$$

or, in terms of $W_T$:

$$\Delta = \partial_{S_0} V_0 = \mathbb{E}(\Phi(\log(S_0)+(r-d-\dfrac{\sigma^2}{2}) T + \sigma W_T) \frac{W_T}{\sigma T S_0} )$$

Let's assume now that $\Phi$ is a forward i.e.:
$$\Phi(x)=e^x$$

Since $\Delta = e^{(r-d)T}$:
$$\mathbb{E}(e^{-\frac{\sigma^2}{2} T + \sigma W_T}\,W_T) = \sigma T$$
\end{section}
\begin{section}{title="Conditional Expectation Method"}
Let's fix a positive number $\tau < T$ and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))$$
\end{section}
\begin{section}{title="Conditional Saltando Expectation Method"}
Let's assume $\tau$ is a stopping time st $\tau < T$ a.s., and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}))$$

Furthemore, if we assume that $X_t$ is a finite activity Levy process, a smart choice of the law of the stopping time can benefit the computation of the inner expectation.
Indeed, if we take $\tau$ as the stopping time corresponding to the last jump of the process, we have that:
$$ X_t | X_{\tau}  \text{ is a Ito process }$$
Why is it useful?\\
Let's assume $X_t$ is a Kou process, then $X_t | X_{\tau}$ is a Brownian motion with a modified drift, hence we can have an analytic formula for the inner expectation.
\end{section}
\begin{section}{title="Vibrato Montecarlo Method"}
Let's fix a positive number $\tau < T$ and let's assume now that the underlying process is a Levy process.\\
Since:
$$V_0=\mathbb{E}(\Phi(X_T))=\mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\mathbb{E}(\mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\Phi(X_T))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_T)|X_{\tau}))=\partial_{\theta} \mathbb{E}(\mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$
$$\partial_{\theta} V_0=\mathbb{E}(\partial_{\theta} \mathbb{E}(\Phi(X_{\tau}+(X_{T}-X_{\tau}))|X_{\tau}))$$

Now we apply the LRM method to the inner expectation, let's call $Y_{T-\tau}= X_{T}-X_{\tau}$, then we have:

$$\partial_{\theta} \mathbb{E}(\Phi(X_{T})|X_{\tau}) = \mathbb{E}(\Phi(X_{T})\partial_{\theta}  \log(f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta))|X_{\tau}) $$

Hence
$$\partial_{\theta} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T})\partial_{\theta}  \log(f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta))|X_{\tau}))$$

or for a generic order n:

$$\partial_{\theta^n} V_0=\mathbb{E}(\mathbb{E}(\Phi(X_{T}) B_n(\{\partial_{\theta^i} \log{f_{X_T-X_{\tau}}(X_T - X_{\tau},\theta)}\}_{i = 1 : n})))$$


As the traditional LR method, this is hard to apply in case the transition density is unknown. This can be easily solved in case of Ito processes where we can approximate the increment as a sum of powers of normal random variables.
In case of more general process like Levy, one needs to invert numerically the density from the characteristic function.

\end{section}

\begin{section}{title="Vibrato Saltando Montecarlo Method"}
Since what we stated above applies also to a generic stopping time $\tau$, we investigate some smart choice of the law of $\tau$ in order to exploit some property of the underlying process.
Let's assume now that $X_t$ is a Finite Activity Levy process and that $\tau$ is the stopping time corresponding to the last jump of the underlying.
We know already that $X_t|X_{\tau}$ is a Ito process.
On such process we are able to apply the standard LRM or VBM since we know the corresponding law of the increments (powers of gaussians).

\end{section}


<!-- ==============================
     GETTING STARTED
     ============================== -->
\begin{section}{title="Getting started"}


In order to get started, just add the package (with **Julia ≥ 1.3**) and

```julia-repl
julia> using PkgPage
julia> newpage()
julia> serve()
```

\\

The `newpage` call will
* generate a `page/` folder in your current directory which contains all the material to generate your site,
    * modify `page/config.md` to change the description, layout, colors etc.,
    * modify `page/index.md` to change the content.
* place a `.github/workflows/DeployPage.yml` to enable the github action that will help deploy your page.

The `serve` call will render your page and make it available for live-preview in your browser when you make modifications to either `config.md` or `index.md`.

\alert{You can specify another folder name via `newpage(path="page2")` but don't forget to modify this in the `DeployPage.yml` as well (in 2 spots).}

\end{section}



<!-- ==============================
     SPECIAL COMMANDS
     ============================== -->
\begin{section}{title="Commands"}

\lead{PkgPage makes a few special commands available to you to simplify the insertion of common useful environments.
}

* [Sections](#com-sections)
* [Figures](#com-figures)
* [Tables](#com-tables)
* [Columns](#com-columns)
* [Maths](#com-maths)
* [Misc](#com-misc)

\label{com-sections}
**Sections**: you can indicate a section as follows:

```plaintext
\begin{section}{opts...}
text
\end{section}
```

where the available options are:

* `title="..."` to add a heading to the section, this should be provided,
* `name="..."` to add a specific name that will appear in the navbar (if not given, it will use the given title).
* `width=8` an integer number controlling the width of the section with respect to the body, increase it to get a wider section.
\label{com-figures}
**Figures**: you can indicate a figure as follows:

```plaintext
\figure{path="...", opts...}
```

where `path` is a valid path to the image (e.g. `/assets/image.png`) and the other available options are:

* `width="..."` to specify the width of the image (e.g. `100%`),
* `alt="..."` to specify the image `alt`,
* `caption="..."` to specify the figure caption if any,
* `style="..."` to add any specific CSS styling to the image (e.g. `border-radius:5px`).

\begin{center}
  \figure{path="/assets/nice_image.jpg", width="100%", style="border-radius:5px;", caption="Panoramic view of the Tara Cathedrals (taken from Wikimedia)."}
\end{center}

\label{com-tables}
**Tables**: you can insert a table as follows:

```plaintext
\begin{table}{opts...}
| Column One | Column Two | Column Three |
|----------: | ---------- |:------------:|
| Row `1`    | Column `2` |              |
| *Row* 2    | **Row** 2  | Column ``3`` |
\end{table}
```

where the available options are:

* `caption="..."`: to add a caption to the table,
* `class="..."`: to add specific bootstrap classes to the table (e.g. `table-striped`),
* `style="..."`: for any further styling.

\begin{center}
\begin{table}{caption="A simple table", class="table-striped"}
| Column One | Column Two | Column Three |
|:---------- | ---------- |:------------:|
| Row `1`    | Column `2` |              |
| *Row* 2    | **Row** 2  | Column ``3`` |
\end{table}
\end{center}

\label{com-columns}
**Columns**: you can declare an environment with columns with:

```plaintext
\begin{columns}
\begin{column}{}
...
\end{column}
\begin{column}{}
...
\end{column}
\end{columns}
```

For instance you can use this to produce:

\begin{columns}
\begin{column}{}
**_Content of a first column here_**

Here is some more content for that first column.
\end{column}
\begin{column}{}
**_Content of a second column here_**

Here is some more content for that second column.
\end{column}
\end{columns}

\\

\label{com-maths}
**Maths**

Just use `$$ ... $$` for display math and  `$ ... $` for inline maths:

$$ w_i = {2 \over (1-x_i^2)[P'_n(x_i)]^2} $$ <!--_-->

where $P_n$ is the $n$-th [Legendre polynomial](https://en.wikipedia.org/wiki/Legendre_polynomials) and $x_i$ is it's $i$-th root.

\label{com-misc}
**Misc**

* you can use `\\` to add some vertical space (skip a line)
\\
* you can use `\begin{center}...\end{center}` to center some content
\begin{center}
_some centered content_
\\\\
\end{center}
* you can use `\style{css}{content}` to get css-styled text for instance `\style{color:red;text-transform:uppercase;}{hello}` gives \style{color:red;text-transform:uppercase;}{hello}
* you can use `\alert{...}` to draw attention to some text via a coloured box.

\alert{this is an alert}

You can also define your own commands which can be as complex as you might want, see the [Franklin docs](https://franklinjl.org) for more information.

\end{section}


<!-- =============================
     SHOWING CODE
    ============================== -->

\begin{section}{title="Showing Code"}

\lead{
    Franklin can run your Julia code on the fly and show the output.
}

**Setting up the environment**: the first step is to ensure that the folder with your source has the proper environment including your package.
To do so, in the Julia REPL, navigate to the source (e.g. `cd("page/")`), activate the environment (e.g. `using Pkg; Pkg.activate()`) and add the package(s) that you need (e.g. `Pkg.add("DataFrames")`).
If you check the status or the Project.toml, you will see that `Franklin` is already in there on top of whatever packages you might have chosen to add.
In our current case:

```
Status `~/.julia/dev/PkgPage/page/Project.toml`
  [a93c6f00] DataFrames v0.21.2
  [713c75ef] Franklin v0.8.2
```

Also add the package in the `DeployPage.yml` e.g. in our case there is:

```julia
Pkg.add(["NodeJS", "DataFrames"]);
```

Once that's set up, you can use "named" code blocks i.e. code blocks that look like

`````
```julia:ex
using DataFrames
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
first(df, 3)
```
`````

where `:ex` is the "named part" (`ex` being the name, which should be unique on the page).

```julia:ex
using DataFrames
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
first(df, 3)
```

You can control the indentation and appearance of the output block in the `config.md` too.

\end{section}


<!-- =============================
     Deploying
    ============================== -->

\begin{section}{title="Deployment"}

\lead{Make your page available online easily by leveraging GitHub Actions and GitHub Pages.}

By following these instructions, the content of the rendered website will be copied to a `gh-pages` branch where it will be deployed by GitHub.
If you would like to deploy the page with your own URL or using something else than GitHub, have a look at the specific instructions further on.

**Adjust DeployPage**: start by checking the `.github/workflows/DeployPage.yml` in particular:
* if you want to use Python or matplotlib, uncomment the relevant lines
* in the `run` block ensure that
    * `NodeJS` and `PkgPage` are added,
    * any packages that your page might rely on are added,
    * the `optimize` call has the appropriate `input` and `output` path (if you're in the default setting, leave as is).

**GitIgnore**: it's important you specify that `page/__site` should be ignored by git and not pushed to your repository otherwise the build process might not work properly. To do so create a file `.gitignore` containing the line

```
page/__site
```

as shown [here](https://github.com/tlienart/PkgPage.jl/blob/cce098535eb95c2c3ba919d605792abfee57710c/.gitignore#L3).

**GitAttributes**: in order for GitHub to ignore `page` folder it the language statistics for your repository, make sure to add a file `.gitattributes` with content

```
page/* linguist-vendored
```

like [this](https://github.com/tlienart/PkgPage.jl/blob/master/.gitattributes).

Now whenever you push changes to the `master` branch of your package, the  build process will be triggered and your page updated and deployed.
**That's it**.

**Avoiding clashes with Documenter.jl**: if you already use [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) you might want your page to be deployed in a specific folder of `gh-pages` as Documenter also generates files in `gh-pages`.

\alert{This will typically not be necessary as the names created by PkgPage and Documenter don't clash, but you might still prefer to avoid mixing the two (in which case, read on).}

you can do so in two steps:

1. change the `run` part of `DeployPage.yml` by specifying the `output` keyword argument  in `PkgPage.optimize` for instance: `PkgPage.optimize(input="page", output="page")`,
1. change the `prepath` in `config.md` to reflect that the base URL will contain that additional folder, for instance `@def prepath = "YourPackage.jl/page"`.

**Use your own URL**: you can usually get host services like Netlify to deploy a specific branch of a GitHub repo, do make sure to set `@def prepath = ""` in your `config.md` though.

If you want to do the deployment without GitHub actions then you will need to:

* ensure you have `purgecss` and `highlights` installed and available to `NodeJS`, the simplest way to do this is to install them via `NodeJS` with

```
using NodeJS;
run(`$(npm_cmd()) install highlight.js`);
run(`$(npm_cmd()) install purgecss`);
```
\\
* run `PkgPage.optimize(input="page", output="")` (adapting `input` as required)
* place the content of `page/__site` wherever your server requires it.

\end{section}
\theoremscripts