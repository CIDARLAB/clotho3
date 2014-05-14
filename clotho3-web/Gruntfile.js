// Generated on 2013-11-18 using generator-angular 0.6.0-rc.1
'use strict';

var shell = require('shelljs');

module.exports = function (grunt) {
  require('load-grunt-tasks')(grunt);
  require('time-grunt')(grunt);

  grunt.initConfig({
    yeoman: {
      // configurable paths
      app: require('./bower.json').appPath || 'app',
      dist: 'dist',
	    api: 'api',
	    command: 'command'
    },
    watch: {
      coffee: {
        files: ['<%= yeoman.app %>/scripts/**/*.coffee'],
        tasks: ['coffee:dist']
      },
      coffeeTest: {
        files: ['test/spec/**/*.coffee'],
        tasks: ['coffee:test']
      },
	    js: {
		    files: ['<%= yeoman.app %>/scripts/**/*.js'],
		    tasks: ['newer:jshint:all'],
		    options: {
			    livereload: true
		    }
	    },
	    jsTest: {
		    files: ['test/spec/**/*.js'],
		    tasks: ['newer:jshint:test', 'karma']
	    },
      compass: {
        files: ['<%= yeoman.app %>/styles/**/*.{scss,sass}'],
        tasks: ['compass:server', 'autoprefixer']
      },
      styles: {
        files: ['<%= yeoman.app %>/styles/**/*.css'],
        tasks: ['copy:styles', 'autoprefixer']
      },
	    gruntfile: {
		    files: ['Gruntfile.js']
	    },
      livereload: {
        options: {
          livereload: '<%= connect.options.livereload %>'
        },
        files: [
          '<%= yeoman.app %>/**/*.html',
          '<%= yeoman.app %>/css/**/*.css',
          '{.tmp,<%= yeoman.app %>}/scripts/**/*.js',
          '<%= yeoman.app %>/images/**/*.{png,jpg,jpeg,gif,webp,svg}'
        ]
      }
    },
	  //todo - update to /css/ not /.tmp/ layout
    autoprefixer: {
      options: ['last 1 version'],
      dist: {
        files: [{
          expand: true,
          cwd: '<%= yeoman.app %>/css',
          src: '**/*.css',
          dest: '<%= yeoman.app %>/css'
        }]
      }
    },
    connect: {
      options: {
        port: 9000,
        // Change this to '0.0.0.0' to access the server from outside.
        hostname: 'localhost',
        livereload: 35729
      },
      livereload: {
        options: {
          open: true,
          base: [
            '.tmp',
            '<%= yeoman.app %>'
          ]
        }
      },
      test: {
        options: {
          port: 9001,
          base: [
            '.tmp',
            'test',
            '<%= yeoman.app %>'
          ]
        }
      },
      dist: {
        options: {
          base: '<%= yeoman.dist %>'
        }
      }
    },
    clean: {
      dist: {
        files: [{
          dot: true,
          src: [
            '.tmp',
	          '<%= yeoman.app %>/css',
            '<%= yeoman.dist %>/*',
            '!<%= yeoman.dist %>/.git*'
          ]
        }]
      },
      server: ['.tmp', 'css']
    },
    jshint: {
      options: {
        jshintrc: '.jshintrc',
        reporter: require('jshint-stylish')
      },
      all: [
        'Gruntfile.js',
        '<%= yeoman.app %>/scripts/**/*.js'
      ],
	    test: {
		    options: {
			    jshintrc: 'test/.jshintrc'
		    },
		    src: ['test/spec/**/*.js']
	    }
    },
    coffee: {
      options: {
        sourceMap: true,
        sourceRoot: ''
      },
      dist: {
        files: [{
          expand: true,
          cwd: '<%= yeoman.app %>/scripts',
          src: '{,*/}*.coffee',
          dest: '.tmp/scripts',
          ext: '.js'
        }]
      },
      test: {
        files: [{
          expand: true,
          cwd: 'test/spec',
          src: '{,*/}*.coffee',
          dest: '.tmp/spec',
          ext: '.js'
        }]
      }
    },

	  // Automatically inject Bower components into the app
	  'bower-install': {
		  app: {
			  html: '<%= yeoman.app %>/index.html',
			  ignorePath: '<%= yeoman.app %>/'
		  }
	  },

    compass: {
      options: {
        sassDir: '<%= yeoman.app %>/styles',
        cssDir: '<%= yeoman.app %>/css',
        generatedImagesDir: '.tmp/images/generated',
        imagesDir: '<%= yeoman.app %>/images',
        javascriptsDir: '<%= yeoman.app %>/scripts',
        fontsDir: '<%= yeoman.app %>/fonts',
        importPath: '<%= yeoman.app %>/bower_components',
        httpImagesPath: '/images',
        httpGeneratedImagesPath: '/images/generated',
        httpFontsPath: '/fonts',
        relativeAssets: false
      },
      dist: {
	      options: {
		      debugInfo : false,
		      generatedImagesDir: '<%= yeoman.dist %>/images/generated'
	      }
      },
      server: {
        options: {
          debugInfo: true
        }
      }
    },
    // not used since Uglify task does concat,
    // but still available if needed
    /*concat: {
      dist: {}
    },*/
	  //fixme - for relative paths of format ../*, revved files don't always work
	  //useminprepare will look in /app/ and usemin will look in /dist/ .. something weird happens. check usemin readme
	  //probably just need to copy over items
    rev: {
      dist: {
        files: {
          src: [
            '<%= yeoman.dist %>/scripts/**/*.js',
            '<%= yeoman.dist %>/styles/**/*.css',
            //'<%= yeoman.dist %>/images/**/*.{png,jpg,jpeg,gif,webp,svg}',
            '<%= yeoman.dist %>/styles/fonts/*'
          ]
        }
      }
    },
    useminPrepare: {
	    html: '<%= yeoman.app %>/index.html', //for single-target only
	    options: {
		    dest: '<%= yeoman.dist %>'
	    }
	    /*dist: {
		    src: '<%= yeoman.app %>/index-dist.html',
	      dest: '<%= yeoman.dist %>'
	    }
	    api: {
		    src: '<%= yeoman.app %>/index-api.html',
		    dest: '<%= yeoman.dist %>'
	    },
	    command: {
		    src: '<%= yeoman.app %>/index-command.html',
		    dest: '<%= yeoman.dist %>'
	    }*/
    },
    usemin: {
      html: ['<%= yeoman.dist %>/**/*.html'],
      css: ['<%= yeoman.dist %>/styles/**/*.css'],
      options: {
        assetsDirs: ['<%= yeoman.dist %>']
      }
    },
    imagemin: {
      dist: {
        files: [{
          expand: true,
          cwd: '<%= yeoman.app %>/images',
          src: '{,*/}*.{png,jpg,jpeg}',
          dest: '<%= yeoman.dist %>/images'
        }]
      }
    },
    svgmin: {
      dist: {
        files: [{
          expand: true,
          cwd: '<%= yeoman.app %>/images',
          src: '{,*/}*.svg',
          dest: '<%= yeoman.dist %>/images'
        }]
      }
    },
    htmlmin: {
      dist: {
	      options: {
		      collapseWhitespace: true,
		      collapseBooleanAttributes: true,
		      removeCommentsFromCDATA: true,
		      removeOptionalTags: true
	      },
        files: [{
          expand: true,
          cwd: '<%= yeoman.dist %>',
          src: ['partials/**/*.html', 'views/**/*.html', 'extensions/**/*.html'],
          dest: '<%= yeoman.dist %>'
        }]
      }
    },
    // Put files not handled in other tasks here
    copy: {
      dist: {
        files: [{
          expand: true,
          dot: true,
          cwd: '<%= yeoman.app %>',
          dest: '<%= yeoman.dist %>',
          src: [
            '*.{ico,png,txt}',
            '.htaccess',
	          '*.html',
            'bower_components/**/*',
            'images/**/*',
            'fonts/*',
	          'extensions/*',
	          'lib/*',
	          'models/*',
	          'partials/**/*',
	          'views/**/*'
          ]
        }, {
          expand: true,
          cwd: '.tmp/images',
          dest: '<%= yeoman.dist %>/images',
          src: [
            'generated/*'
          ]
        }]
      },
      styles: {
        expand: true,
        cwd: '<%= yeoman.app %>/styles',
        dest: 'css/',
        src: '**/*.css'
      }
    },
    concurrent: {
      server: [
        'coffee:dist',
        'compass:server',
        'copy:styles'
      ],
      test: [
        'coffee',
        'compass',
        'copy:styles'
      ],
      dist: [
        'coffee',
        'compass:dist',
        'copy:styles',
        //'imagemin',
        'svgmin',
        'htmlmin'
      ]
    },
    karma: {
      unit: {
        configFile: 'karma.conf.js',
        singleRun: true
      },
	    api : {
		    configFile: 'karma-api.conf.js',
		    singleRun: true
	    },
	    command : {
		    configFile: 'karma-command.conf.js',
		    singleRun: true
	    }
    },
    cdnify: {
      dist: {
        html: ['<%= yeoman.dist %>/*.html']
      }
    },
    ngmin: {
      dist: {
        files: [{
          expand: true,
          cwd: '.tmp/concat/scripts',
          src: '*.js',
          dest: '.tmp/concat/scripts'
        }]
      }
    },
    uglify: {
      dist: {
        /*
        //this was the default, but we broke up the scripts into multiple modules... will still use generated block
        files: {
          '<%= yeoman.dist %>/scripts/scripts.js': [
            '<%= yeoman.dist %>/scripts/scripts.js'
          ]
        }
        */
      }
    },
	  //grunt-conventional-changelog
	  //todo - use
	  changelog: {
		  options: {
			  // Task-specific options go here.
		  }
	  },
	  shell: {
		  mongo: {
			  command: "sh startMongoIfNotRunning.sh",
			  options: {
				  async: true
			  }
		  },
		  clothoCleanBuild: {
			  command: 'cd ..; mvn clean install'
		  },

		  /*
		   * The following commands for running the server require mvn command line tools installed:
	     * export PATH=/usr/local/apache-maven-3.1.1/bin:$PATH
			 */

		  /**
		   * Run the test server, from app (client source files)
		   */
		  clothoTestServer: {
			  /*

			  */
			  command: 'cd ..; mvn "-Dexec.args=-Dloglevel="OFF" -classpath %classpath org.clothocad.core.util.ClothoTestEnvironment -clientdirectory clotho3-web/app" -Dexec.executable=java -Dexec.classpathScope=test process-classes org.codehaus.mojo:exec-maven-plugin:1.2.1:exec',
			  options: {
				  async: true
			  }
		  },
			clothoTestServerVerbose: {
			  command: 'cd ..; mvn "-Dexec.args= -classpath %classpath org.clothocad.core.util.ClothoTestEnvironment -clientdirectory clotho3-web/app" -Dexec.executable=java -Dexec.classpathScope=test process-classes org.codehaus.mojo:exec-maven-plugin:1.2.1:exec',
			  options: {
				  async: true
			  }
		  },
		  /**
		   * Run the Production server, from app (client source files)
		   */
		  clothoAuthoringAppServer: {
			  command: 'cd ..; mvn "-Dexec.args=-Dloglevel="OFF" -classpath %classpath org.clothocad.core.util.ClothoAuthoringEnvironment -clientdirectory clotho3-web/app" -Dexec.executable=java -Dexec.classpathScope=test process-classes org.codehaus.mojo:exec-maven-plugin:1.2.1:exec',
			  options: {
				  async: true
			  }
		  },
		  /**
		   * Run the Production server, from dist (client compiled files)
		   */
		  clothoProdServer: {
			  command: 'cd ..; mvn "-Dexec.args=-Dloglevel="OFF" -classpath %classpath org.clothocad.core.util.ClothoAuthoringEnvironment -clientdirectory clotho3-web/dist" -Dexec.executable=java -Dexec.classpathScope=test process-classes org.codehaus.mojo:exec-maven-plugin:1.2.1:exec',
			  options: {
				  async: true
			  }
		  }
	  },
	  open : {
		  dev : {
			  path: 'https://localhost:8443/',
			  app: 'Google Chrome'
		  }
	  },
	  //todo - use in workflow
	  "merge-conflict": {
		  files: [
			  '**/*'
		  ]
	  },
	  // Removes unused CSS
	  //future - integrate when handle dynamically added css
	  uncss: {
		  dist: {
			  files : {
				  '<%= yeoman.dist %>/styles/main.css': ['<%= yeoman.dist %>/*.html', '<%= yeoman.dist %>/views/**/*','<%= yeoman.dist %>/partials/**/*']
			  }
		  }
	  },
		//note - this is for versioning
	  bump: {
		  options: {
			  files: ['package.json'],
			  commit: true,
			  createTag: true,
			  push: true,
			  pushTo: 'origin'
		  }
	  },
	  processhtml : {
		  options : {
			  commentMarker : 'process', //don't want to use default 'build' bc usemin,
			  stripUnparsed : true
		  },
		  api : {
			  files : {
				  '<%= yeoman.dist %>/index-api.html': ['<%= yeoman.dist %>/index.html']
			  }
		  },
		  command : {
			  files : {
				  '<%= yeoman.dist %>/index-command.html': ['<%= yeoman.dist %>/index.html']
			  }
		  },
		  dist : {
			  files : {
				  '<%= yeoman.dist %>/index.html': ['<%= yeoman.dist %>/index.html']
			  }
		  }
	  }
  });


	grunt.registerTask('update:npm', 'Update NPM packages.', function () {
		shell.exec('npm install');
	});


	grunt.registerTask('serve', function (target) {
    if (target === 'dist') {
      return grunt.task.run(['build', 'connect:dist:keepalive']);
    }

    grunt.task.run([
      'clean:server',
      'concurrent:server',
      'autoprefixer',
      'connect:livereload',
      'watch'
    ]);
  });

	grunt.registerTask('dev', [
		'shell:mongo',
		'shell:clothoTestServer',
		'clean:server',
		'concurrent:server',
		'autoprefixer',
		'open',
		'watch'
	]);

	grunt.registerTask('authoring', [
		'shell:mongo',
		'shell:clothoAuthoringAppServer',
		'clean:server',
		'concurrent:server',
		'autoprefixer',
		'open',
		'watch'
	]);

	grunt.registerTask('prod', [
		'shell:mongo',
		'shell:clothoProdServer',
		'clean:server',
		'concurrent:server',
		'autoprefixer',
		'open',
		'watch'
	]);

  grunt.registerTask('test', [
    'clean:server',
    'concurrent:test',
    'autoprefixer',
    'connect:test',
    'karma'
  ]);

  grunt.registerTask('build', [
    'clean:dist',
    'useminPrepare',
    'concurrent:dist',
    'autoprefixer',
    'concat',
    'ngmin',
    'copy:dist',
    //'cdnify',
    'cssmin',
    'uglify',
    //'rev',
    'usemin',
	  'processhtml',
	  'htmlmin'
  ]);

	grunt.registerTask('build-uncss', [
		'build',
		'uncss'
	]);

  grunt.registerTask('default', [
    'newer:jshint',
    'test',
    'build'
  ]);
};
