import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("imu_gnss_data.csv")  # Replace with actual file name

# Convert time to relative (starting from zero)
df["imu_time"] = df["imu_time"] - df["imu_time"].min()

# Plot gyroscope
plt.figure(figsize=(10, 4))
plt.plot(df["imu_time"], df["gyrox"], label="Gyro X")
plt.plot(df["imu_time"], df["gyroy"], label="Gyro Y")
plt.plot(df["imu_time"], df["gyroz"], label="Gyro Z")
plt.title("Gyroscope over Time")
plt.xlabel("Time [s]")
plt.ylabel("Gyro [rad/s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot magnetometer
plt.figure(figsize=(10, 4))
plt.plot(df["imu_time"], df["magx"], label="Mag X")
plt.plot(df["imu_time"], df["magy"], label="Mag Y")
plt.plot(df["imu_time"], df["magz"], label="Mag Z")
plt.title("Magnetometer over Time")
plt.xlabel("Time [s]")
plt.ylabel("Magnetic Field [uT?]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

