# DiSCo: Distributional Synthetic Controls for Stata <img src="figures/logo.png" align="right" alt="" width="155" />

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/Davidvandijcke/DiSCos_stata)](https://github.com/Davidvandijcke/DiSCos_stata/releases)
[![GitHub last commit](https://img.shields.io/github/last-commit/Davidvandijcke/DiSCos_stata.svg)](https://github.com/Davidvandijcke/DiSCos_stata/commits/main)
![GitHub issues](https://img.shields.io/github/issues/Davidvandijcke/DiSCos_stata)

🕺 Let's get DiSCo-ing with distributions! This package implements Distributional Synthetic Controls (DiSCo) in Stata, following the methodology developed in Gunsilius (2023). While traditional synthetic controls focus on matching means, DiSCo takes it to the next level by matching entire distributions.

## 🚀 Features

- 📊 Match entire outcome distributions, not just means
- 🎯 Estimate distributional treatment effects
- 📈 Generate ~~beautiful~~ visualizations
- 🔄 Bootstrap inference and permutation tests
- 🎨 Flexible estimation approaches (quantile-based or CDF-based)

## 💿 Installation

The package will soon be available on SSC. For now, you can install directly from GitHub:

```stata
net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/main/src/") replace
```

For the development version:
```stata
net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/dev/src/") replace
```

## 🎵 Quick Start

```stata
* Load your data
use mydata.dta, clear

* Run DiSCo
disco outcome unit time, idtarget(1) t0(10) graph 

* Get summary statistics
disco_estat
```

## 📚 Documentation

After installation, access the help files in Stata:
```stata
help disco
help disco_estat
help disco_plot
```

## 🎪 In Development

This package is actively being developed! Future releases may include (depending on my mood and the weather):
- Additional visualization options
- More flexible weighting schemes
- [Multivariate distributional synthetic controls... (scary)](https://www.jmlr.org/papers/v25/23-0708.html)

Comments, suggestions, and bug reports are very welcome! Please open an issue on GitHub or email dvdijcke@umich.edu.

## 🎵 Related Packages

Also check out the R version of DiSCo: [![CRAN](https://www.r-pkg.org/badges/version/DiSCos)](https://cran.r-project.org/package=DiSCos)

## 📝 Citation

If you use this package, please cite:

```bibtex
@article{gunsilius2023distributional,
  title={Distributional Synthetic Controls},
  author={Gunsilius, Florian F},
  journal={Econometrica},
  volume={91},
  number={3},
  pages={1105--1117},
  year={2023}
}
@article{van2024return,
  title={Return to Office and the Tenure Distribution},
  author={Van Dijcke, David and Gunsilius, Florian and Wright, Austin L},
  journal={University of Chicago, Becker Friedman Institute for Economics Working Paper},
  number={2024-56},
  year={2024}
}

```

## 💃 Why do the DiSCo?

Because averages are so yesterday! If you're feeling like doing the robot (you're thinkin' bout that synthetic control), and you want your treatment effects to be as varied as your moves on the dance floor, you'll need DiSCo to capture all that distributional groove. 

## 🎉 License

This project is licensed under the MIT License.

---
Made with 🕺 and 💃 by [David Van Dijcke](https://www.davidvandijcke.com/), with a stamp of approval from [Florian Gunsilius](https://www.floriangunsilius.com/)