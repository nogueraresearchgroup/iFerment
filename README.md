# iChainElongateSugars

**iChainElongateSugars** is a single-cell metabolic model for a general bacterium capable of sugar chain-elongation. It is an amended version of Matthew Scarborough’s unpublished iFerment181 model, which includes biomass growth, for which all vestiges remain, but may be disabled via commenting out the code.

### Dependencies:

- Anaconda distribution of Python 3.
- COBRApy
- Various other libraries, see `environment.yml` for a complete list.

### Operation:

1. If you do not have an Anaconda installation, install Anaconda by following the guide [here](https://docs.anaconda.com/anaconda/user-guide/getting-started/).
2. Clone this repository.
   ```bash
   git clone <INSERT URL HERE>
   ```
3. Create a new Anaconda environment, installing the needed libraries:
   ```bash
   conda env create --file environment.yml
   ```
4. Run the script from the terminal.
   ```bash
   python3 iChainElongateSugars.py
   ```

## Model Organization

### Part 1: Reactor Parameters:

The values in this section are only relevant if “Transenergetics == True” a.k.a. Transport Energetics, i.e., if we want to impose our metabolomics data onto the model, we can manually do this in this section. If “Transenergetics == False” and we do not impose such data, we get different results. The next time these values show up is in line 3518 (bottom of Part 2) where the *imposing* part happens.  

- Experimental extracellular steady-state values
- Assumed intracellular steady-state values
- Textbook thermodynamics used in membrane transport



### Part 2: Build the Model:

> NOTE: Reactions in CobraPy do not have to be balanced. In this section, there are manual additions of metabolites, reactions, and their definitions in accordance with CobraPy rules (this is the bulk of the model).
>
> “Exchange” reactions account for metabolites that can be extracellular. 
>
> “Import/export” reactions account for metabolites that can leave the system provided their Exchange results in some flux and theoretically reenter the system (the reentering part is not necessarily true because of the default bounds set in the model). 

- Carbon metabolism
- Adjacent metabolism
- Reactions for biomass
- Import/export reactions  
- ATP Transport reactions – these are necessary for “Transenergetics“ constraints
- Transport energy – this is where the flux of ATP hydrolysis attributed to each end metabolite is set to the following:  
  - membrane 	transport thermodynamic-based coefficient*(flux attributed to 	export of metabolites)

### Part 3: Substrate Uptake:

Here, you are effectively declaring the specific maximum utilization rate of any metabolite for which you’ve provided an Exchange reaction – you are setting an upper bound on the flux (mmol/time/1 g DCW) of this reaction. This means reported fluxes will be in the same units . . . which makes interpreting the BIOMASS flux tricky – I’ll have to run some thoughts by you.

### Part 4: Set Additional Constraints:

Here’s a list of reactions one can “knock_out” to see what happens. This is also where the model’s objective is set, where several documentative values are calculated, and where the generation of output reports are commenced.

- Reaction knockouts
- Thermo. reports
- Print outputs

---

Amendment 1: Essential Proteins: This is a for-loop that goes deletes each reaction, then prints the ATP production after deletion of that reaction. You can see some of the genes in here are essential for ATP production – I will explore the substrate dependence of this. 