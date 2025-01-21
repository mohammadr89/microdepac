extern "C" {
    void calc_lai_wrapper(
        int day_of_year,
        float lat,
        int lu,
        float* lai,
        float* sai,
        float* laimax);
}
