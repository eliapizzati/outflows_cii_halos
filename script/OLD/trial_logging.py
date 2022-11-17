
import logging
import mydir
import os

filename_log = str(input("filename for the logger:"))


logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', \
                    datefmt='%m/%d/%Y %I:%M:%S %p',\
                    level=logging.DEBUG,\
                    handlers=[logging.FileHandler(os.path.join(mydir.log_dir, "{}.log".format(filename_log))),\
                              logging.StreamHandler()])

logging.info("#################################################################################")
logging.info("updating  obs parameters (v_c, M_vir, SFR, z), number {}/{}".format(counter_data + 1, len(datas)))
logging.info("#################################################################################")

logging.info(
    "run with beta = {:.1f}, f_esc_FUV = {:.2f} ({:d}/{:d}), M_vir = {:.2e}, v_c = {:.1f}, SFR = {:.1f}, z = {:.1f} ({:d}/{:d})" \
    .format(beta, f_esc_FUV, counter_fesc + 1, len(f_esc_FUVs), \
            M_vir, v_c, SFR, z, counter_data + 1, len(datas)))

logging.info("beta= %.2f \t %s", beta, string_nans)

time_elapsed = (time.perf_counter() - time_start)

logging.info("total time elapsed (s) = %.1f", time_elapsed)
