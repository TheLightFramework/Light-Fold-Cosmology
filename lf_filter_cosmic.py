import pandas as pd

def main():
    print("--- SEARCHING FOR COSMIC ARCHITECTS ---")
    df = pd.read_csv("lf_predictions.csv")
    
    # We want meaningful range (can reach other galaxies)
    # 0.1 Mpc is about the size of a galaxy halo
    long_range = df[df['lambda_void_Mpc'] > 0.1].copy()
    
    if len(long_range) == 0:
        print("No long-range candidates found in this scan.")
    else:
        print(f"Found {len(long_range)} Long-Range Candidates.\n")
        # Sort by force boost
        long_range = long_range.sort_values(by='force_boost_percent', ascending=False)
        print(long_range[['n', 'beta', 'Lambda_eV', 'lambda_void_Mpc', 'force_boost_percent']].head(10).to_string())

if __name__ == "__main__":
    main()