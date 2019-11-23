<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App;

use App\Exceptions\CommandException;
use App\Exceptions\ProcessingJobException;
use Symfony\Component\Process\Exception\ProcessFailedException;
use Symfony\Component\Process\Process;

final class Utils
{

    /**
     * Runs a shell command and checks for successful completion of execution
     *
     * @param array         $command
     * @param string|null   $cwd
     * @param int|null      $timeout
     * @param callable|null $callback
     *
     * @return string|null
     */
    public static function runCommand(
        array $command,
        ?string $cwd = null,
        ?int $timeout = null,
        ?callable $callback = null
    ): ?string {
        $process = new Process($command, $cwd, null, null, $timeout);
        $process->run($callback);
        if (!$process->isSuccessful()) {
            throw new ProcessFailedException($process);
        }

        return $process->getOutput();
    }

    /**
     * Map command exception to message
     *
     * @param \Symfony\Component\Process\Exception\ProcessFailedException $e
     * @param array                                                       $errorCodeMap
     *
     * @return \App\Exceptions\ProcessingJobException
     */
    public static function mapCommandException(
        ProcessFailedException $e,
        array $errorCodeMap = []
    ): ProcessingJobException {
        $code = $e->getProcess()->getExitCode();
        if (isset($errorCodeMap[$code])) {
            return new ProcessingJobException($errorCodeMap[$code], $code, $e);
        }

        return new ProcessingJobException($e->getMessage(), $code, $e);
    }

}
