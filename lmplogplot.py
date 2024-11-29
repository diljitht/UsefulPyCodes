import sys
import re
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

matplotlib.rcParams.update(
    {
        "font.size": 14,
        "font.family": "STIX",
        "font.weight": "bold",
        "mathtext.fontset": "stix",
    }
)


class Customlabel:
    def __init__(
        self, ax, axis="y", label="", unit="", color="k", tickwidth=2, axeswidth=2
    ):
        self.ax = ax
        self.axis = {"y": ax.yaxis, "x": ax.xaxis}[axis]
        self.label = label
        self.unit = unit
        self.color = color
        self.tickwidth = tickwidth
        self.axeswidth = axeswidth
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        self.axis.set_major_formatter(formatter)
        self.axis.set_minor_locator(ticker.AutoMinorLocator())
        self.axis.label.set_color(self.color)
        [sp.set_linewidth(self.axeswidth) for sp in self.ax.spines.values()]
        self.axis.set_tick_params(
            which="both", width=self.tickwidth, direction="in", colors=self.color
        )
        self.ax.callbacks.connect(axis + "lim_changed", self.update)
        self.ax.figure.canvas.draw()
        self.update(None)

    def update(self, lim):
        fmt = self.axis.get_major_formatter()
        self.axis.offsetText.set_visible(False)
        self.axis.set_label_text(
            self.label + " [" + fmt.get_offset() + self.unit + "]")


# Function to parse LAMMPS log file and return data for each run
def parse_lammps_log(filename):
    runs_data = []
    current_data = {}
    in_thermo = False

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace

            # Detect thermo output based on "Step" header
            if "Step" in line:
                if (
                    current_data
                ):  # If there's existing data, save it before starting a new run
                    runs_data.append(current_data)
                in_thermo = True
                headers = line.split()  # Read headers into a list
                current_data = {
                    header: [] for header in headers
                }  # Initialize new data structure
                continue

            if in_thermo:
                # Read the values in the thermo data
                values = line.split()
                if len(values) > 0:  # Ensure there's at least one value (Step)
                    try:
                        # First column is an integer (Step)
                        step = int(values[0])
                    except ValueError:
                        continue  # Skip this line if it's not a valid integer

                    current_data["Step"].append(step)

                    # Append other values based on headers
                    for i in range(1, len(values)):
                        header = headers[i]
                        try:
                            value = float(values[i])
                            current_data[header].append(value)
                        except ValueError:
                            continue  # Skip invalid values

    if current_data:  # Save the last run data
        runs_data.append(current_data)

    return runs_data


# Function to determine the unit based on regex matching of key names
def get_unit_for_key(key):
    # List of regex patterns and corresponding units
    patterns_to_units = [
        (r"mass", r"$g/mol$"),  # Mass in gram/mole
        (r"distance", r"${\AA}$"),  # Distance in Angstrom
        (r"time", r"$fs$"),  # Time in femtoseconds
        (r"energy|eng|toteng", r"$kcal/mol$"),  # Energy in kilo-calorie/mole
        (r"temp", r"$K$"),  # Temperature in Kelvin
        (r"press", r"$atm$"),  # Pressure in atmospheres
        (r"density", r"$g/cm^{3}$"),  # Density in g/cm^3
        (r"volume", r"$\AA^{3}$"),  # Volume in cubic angstrom
        (r"force", r"$(kcal/mol)/\AA$"),
        (r"msd", r"$\AA^{2}/fs$"),
        (r"Rg", r"$\AA^{2}$"),
        # Add more patterns and units as needed
    ]

    # Iterate over the patterns and return the matching unit
    for pattern, unit in patterns_to_units:
        if re.search(pattern, key, re.IGNORECASE):  # Case-insensitive matching
            return unit

    return "unknown"  # Default if no match


# Function to display available quantities and let the user select them
def user_select_quantities(available_keys):
    print("\nAvailable quantities to plot:")
    for i, key in enumerate(available_keys):
        unit = get_unit_for_key(key)  # Get unit using regex matching
        print(f"{i + 1}. {key} ({unit})")

    # Ask the user to input the numbers corresponding to their selection
    selected_numbers = input(
        "\nEnter the numbers of the quantities to plot (comma-separated): "
    )
    selected_numbers = selected_numbers.split(",")

    # Convert the selected numbers to a list of keys
    selected_keys = [
        available_keys[int(num.strip()) - 1]
        for num in selected_numbers
        if num.strip().isdigit()
    ]

    return selected_keys


# Function to plot data from the log file
def plot_data(runs_data):
    colorset = [
        "red",
        "blue",
        "green",
        "purple",
        "orange",
        "cyan",
        "magenta",
        "yellow",
        "black",
        "brown",
    ]

    for run_index, data in enumerate(runs_data):
        available_keys = [
            key for key in data.keys() if key != "Step" and len(data[key]) > 0
        ]

        if not available_keys:
            print(f"No data to plot for run {run_index + 1}.")
            continue

        # Let the user select which quantities to plot
        selected_keys = user_select_quantities(available_keys)

        if not selected_keys:
            print(f"No valid quantities selected for run {run_index + 1}.")
            continue

        # Create subplots based on selected quantities
        num_plots = len(selected_keys)
        fig, axs = plt.subplots(num_plots, 1, figsize=(15, 5 * num_plots))

        # Handle the case where there's only one subplot (axs is not an array)
        if num_plots == 1:
            axs = [axs]

        # Plot each selected quantity against Time in different subplots
        for idx, (ax, key) in enumerate(zip(axs, selected_keys)):
            color = colorset[
                idx % len(colorset)
            ]  # Cycle through the custom list of colors
            cumulative_avg = np.cumsum(
                data[key]) / np.arange(1, len(data[key]) + 1)
            ax.plot(
                np.array(range(len(data["Step"]))[:: len(data[key]) // 100]) * dt,
                data[key][:: len(data[key]) // 100],
                label=key,
                color=color,
                linewidth=2,
                # linestyle='None',
                linestyle="--",
                marker="o",
            )
            # ax.plot(
            #     np.array(data["Step"]) * dt,
            #     cumulative_avg,
            #     label="cumulative_avg of " + key,
            #     color=color,
            #     linewidth=2,
            #     linestyle='-',
            # )
            ax.grid()
            # Set the x-axis label only for the last subplot
            if idx == num_plots - 1:
                custom_xlabel = Customlabel(
                    ax,
                    label=r"Time",
                    unit=r"$fs$",
                    axis="x",
                    color="black",
                    tickwidth=2,
                    axeswidth=2,
                )
            else:
                # Hide only the labels of the x-ticks but keep the ticks and minor ticks
                ax.tick_params(
                    which="both", width=2, direction="in", labelbottom=False
                )  # This keeps the ticks but hides the labels
                ax.minorticks_on()  # Ensure that minor ticks are enabled for all subplots
            custom_ylabel = Customlabel(
                ax,
                label=f"{key}",
                unit=f"{get_unit_for_key(key)}",
                axis="y",
                color="black",
                tickwidth=2,
                axeswidth=2,
            )

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python " + sys.argv[0] + " <input_file.log>")
    else:
      # Parse the LAMMPS log file
      try:
        log_filename = sys.argv[1]
      except FileNotFoundError:
        print(f"Error: The file '{log_filename}' was not found.")
        exit()

    dt = 1.0
    runs_data = parse_lammps_log(log_filename)
    
    # Plot the data
    plot_data(runs_data)
