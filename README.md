# Feature Characterization
The algoritm is based on the definitions in ISO 21920-2.
With appropriate permissions, the iso can be viewed here: https://nautos.de/3PI/search/item-detail/DE30091970

A documentation of the algorihtm can be found here: https://www.overleaf.com/4687699647jtrjxjypwqmz#0ccd85

## Watershed Segmentation
In the following are some illustrations to show idea of watershed segemenation. For more information see documentation.

### Method to determine watersheds in 2.5D data set
<div align="center">
<video controls src="data/figures for readme/animation.mp4"></video>
</div>

### Method transferred to 2D data set
<div align="center">
<img width="720" src="data/figures for readme/animation.gif" />
</div>

## Usage of feature characterization
The Convention is summarized in the following figure:
<div align="center">
<img width="720" src="data/figures for readme/FC_Convention.png" />
</div>
The functionality can be tested directly using minimal_example.m with or "minimal_example.py" an editable dummy profile. Alternatively, there is a GUI for Matlab "GUI.mlapp" where, for example, the profiles from "data/profles for case studies" can be loaded and the algorithm applied by varying the various input arguments.
<div align="center">
<img width="720" src="data/figures for readme/GUI.png" />
</div>

## Preliminaries Matlab
Add "FC_Functions"-folder to search path of Matlab
```
addpath(*path to FC_Functions*)
```
to permanently save the path
```
save path
```
##### Version
- MATLAB 2017a and higher

## Preliminaries python

Prepare virtual environment, install package and requirements for `matplotlib` \& `jupyter notebook`.

**Windows**
```bash
python -m venv venv
venv\Scripts\activate
pip install python\.
pip install -r python\requirements.txt
```

**Unix**
```bash
python -m venv venv
source venv/bin/activate
pip install python/.
pip install -r python/requirements.txt
```

##### Version
python>=3.11