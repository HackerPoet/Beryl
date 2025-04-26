const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const exe_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
    });

    const exe = b.addExecutable(.{
        .name = "Beryl",
        .root_module = exe_mod,
    });

    const eigen = b.dependency("Eigen", .{
        .target = target,
        .optimize = optimize,
    });

    exe.addCSourceFiles(.{
        .root = b.path(""),
        .files = &.{
        "cmdLine.cpp",
        "globals.cpp",
        "main.cpp",
        "permutation.cpp",
        "solver.cpp",
        "symmetry.cpp",
        "util.cpp",
        },
        .flags = &.{"-std=c++17"},
    });

    exe.addIncludePath(eigen.path(""));

    exe.linkLibCpp();

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);
}
