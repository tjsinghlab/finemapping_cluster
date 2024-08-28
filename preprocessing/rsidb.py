import sqlite3
import yaml
import timeit

class RSIdb:
    def __init__(self, db_path):
        try:
            self.conn = sqlite3.connect(db_path)
            self.cursor = self.conn.cursor()
        except sqlite3.Error as e:
            raise ConnectionError(f"Failed to connect to the database: {e}")

    def query(self, query, params=()):
        if not query.strip():
            raise ValueError("Query string cannot be empty.")
        
        try:
            self.cursor.execute(query, params)
            return self.cursor.fetchall()
        except sqlite3.OperationalError as e:
            raise ValueError(f"Database operation failed: {e}")
        except Exception as e:
            raise e
    
    def lift_rsid_range(self, rsid_low, rsid_high, build="hg19", **kwargs):
        if build not in ["hg19", "hg38"]:
            raise ValueError("Invalid build. Please choose either 'hg19' or 'hg38'.")

        # Get the integer in the rsid, rs0001234 -> 1234
        try:
            rsid_low = int(rsid_low[2:])
            rsid_high = int(rsid_high[2:])
        except ValueError:
            raise ValueError("Invalid rsid format. Please ensure it starts with 'rs' followed by digits.")

        if rsid_low > rsid_high:
            raise ValueError("Invalid rsid range: rsid_low must be less than or equal to rsid_high.")

        query_str = f"""
            SELECT ID, CHROM, POS, REF, ALT
            FROM {build}
            WHERE rsid >= ? AND rsid <= ?;
        """

        verbose = kwargs.get('verbose', False)
        if verbose:
            print("Executing query:")
            print(query_str)
            print("With parameters:", (rsid_low, rsid_high))

        return self.query(query_str, (rsid_low, rsid_high))

    def lift_rsid(self, rsid, build="hg19", **kwargs):
        if build not in ["hg19", "hg38"]:
            raise ValueError("Invalid build. Please choose either 'hg19' or 'hg38'.")

        # Get the integer in the rsid, rs0001234 -> 1234
        try:
            rsid = int(rsid[2:])
        except ValueError:
            raise ValueError("Invalid rsid format. Please ensure it starts with 'rs' followed by digits.")

        query_str = f"""
            SELECT ID, CHROM, POS, REF, ALT
            FROM {build}
            WHERE rsid == ?;
        """

        verbose = kwargs.get('verbose', False)
        if verbose:
            print("Executing query:")
            print(query_str)
            print("With parameters:", (rsid,))

        return self.query(query_str, (rsid,))

    def close(self):
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()

if __name__ == "__main__":
    with open("config.yaml", "r") as file:
        config = yaml.safe_load(file)

    db_path = config['powerlift']['db_path']
    db = RSIdb(db_path)

    # def benchmark_function():
    #     for _ in range(10000):
    #         db.lift_rsid(rsid="rs982", build="hg19", verbose=False)

    # # Run the benchmark
    # number_of_runs = 3
    # execution_time = timeit.timeit(benchmark_function, number=number_of_runs)

    # print(f"Total execution time for {number_of_runs} run(s) of 10,000 queries: {execution_time:.4f} seconds")
    # print(f"Average time per run: {execution_time/number_of_runs:.4f} seconds")
    # print(f"Average queries per second: {10000/(execution_time/number_of_runs):.1f}")

    # open file and store the data as a list
    with open("/gpfs/commons/home/sfriedman/projects/anne_liftover/zeng_rsids_GRCh38.tsv", "r") as file, open("/gpfs/commons/home/sfriedman/projects/anne_liftover/zeng_rsids_GRCh38_expanded.tsv", "w") as output:
        next(file)
        rsids = file.readlines()
        rsids = [rsid.strip() for rsid in rsids]
        print(rsids)
        output.write("rsid\tchromosome\tpos\tpos\n")
        for rsid in rsids:
            res = db.lift_rsid(rsid, build="hg38", verbose=False)
            if res:
                sid, chromosome, pos, ref, alt = res[0]
                output.write(f"{chromosome}:{pos}-{pos}\n")
            else:
                print(f"rsid {rsid} not found")
    

    db.close()