#################################################
##                                             ##
##  Aerodynamic Simulation Automation Program  ##
##                                             ##
#################################################

#################### Overview ###########################

This Python program aims to automate the process of aerodynamic simulation using the OpenFOAM software. It provides a user-friendly interface for generating and refining mesh grids, guiding users through simulation steps from pre-processing to post-processing. By leveraging Python's libraries and a clear structure, the program enhances accessibility to aerodynamic simulation techniques, making it suitable for both beginners and experts in the field.

##################### Requirements ######################

To run this program, you need to have the following installed:

- OpenFOAM 2.2.12: Ensure that you have OpenFOAM version 2.2.12 installed on your system. If you're using a different version, you'll need to modify the version in the standard/incompressible/run_simulation directory to match the version you're using.

- Python: The program is written in Python, so you need to have Python installed on your system. Python 3.0 is recommended.

- tkinter: This program's graphical interface is built using the tkinter module, which is included in standard Python installations.

numpy: Install numpy using the following command:
pip install numpy

matplotlib: Install matplotlib using the following command:
pip install matplotlib

pandas: Install pandas using the following command:
pip install pandas

openpyxl: Install openpyxl using the following command:
pip install openpyxl

tqdm: Install tqdm using the following command:
pip install tqdm


##################### Installation #####################

Install OpenFOAM: Follow the installation instructions for OpenFOAM 2.2.12 specific to your operating system. You can download OpenFOAM 2.2.12 from openfoam.org (https://www.openfoam.com/news/main-news/openfoam-v2212)

Clone Repository: Clone this repository to your local machine using the following command:

bash

git clone https://github.com/your-username/aerodynamic-simulation.git
Replace your-username with your GitHub username.

Usage
Navigate to Directory: Open a terminal or command prompt and navigate to the directory where you cloned the repository.

Run the Program: Run the Python script aero_sim.py using the following command:


python Run.py


##################### Follow Instructions #####################

The program will prompt you with instructions for each step of the simulation process. Follow the on-screen prompts to generate and refine mesh grids, perform simulations, and analyze results.

Select Coordinate Document or Specify NACA Profile: After launching the program, you will be prompted to either select a document containing coordinates (default is 100 points) or specify a desired NACA profile with 4 digits.

Mesh Refinement: While the default mesh can be used for simulation, it is recommended to refine the mesh based on the complexity of the simulation, such as high velocities, compressible behavior, etc.

Provide Simulation Properties: Input simulation properties such as velocity, temperature, pressure, etc.

**Note:**  When specifying multiple angles, the program may take longer to complete due to increased computational complexity.

**Note:** To alter the number of iterations, navigate to one of the standard models in the system folder, open the controlDict document, and modify endTime or deltaT.

**Note:** Os resultados das simulações serão salvos na pasta "resultados", e os gráficos serão salvos na pasta "graficos".


##################### Contributions #####################

Contributions to this project are welcome. If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request on the GitHub repository.

License
This program is licensed under the MIT License. See the LICENSE file for details.

Contact

For any questions or inquiries, feel free to contact the project maintainer at 180078330@aluno.unb.br.