<?php

use Illuminate\Database\Migrations\Migration;
use Illuminate\Database\Schema\Blueprint;
use Illuminate\Support\Facades\Schema;

class CreateJobsTable extends Migration
{
    /**
     * Run the migrations.
     *
     * @return void
     */
    public function up()
    {
        Schema::create('jobs', function (Blueprint $table) {
            $table->bigIncrements('id');
            $table->string('job_type');
            $table->enum('status', ['queued', 'processing', 'completed', 'failed']);
            $table->json('job_input');
            $table->json('job_output');
            $table->text('log');
            $table->unsignedBigInteger('user_id')->index();
            $table->foreign('user_id', 'user_id_to_user_foreign_key')->references('id')->on('users')
                  ->onDelete('cascade')->onUpdate('cascade');
            $table->timestamps();
        });
    }

    /**
     * Reverse the migrations.
     *
     * @return void
     */
    public function down()
    {
        Schema::dropIfExists('jobs');
    }
}
