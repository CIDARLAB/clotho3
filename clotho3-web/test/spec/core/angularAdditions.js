describe('clotho.core angularAdditions', function() {
	var $rootScope;

	beforeEach(function() {
		module('clotho.angularAdditions');

		inject(function (_$rootScope_) {
			$rootScope = _$rootScope_;
		});
	});

	describe('($rootScope.$safeApply)', function () {
		it('should be defined', function () {
			expect(angular.isFunction($rootScope.$safeApply)).toBe(true);
		});

		it('should force a $digest', function () {
			//todo
		});

		it('should work in a $digest', function () {
			//todo
		});
	});


	describe('(angular global)', function () {

		describe('#angular.isEmpty', function () {
			it('should be defined', function () {
				expect(angular.isFunction(angular.isEmpty)).toBe(true);
			});

			it('should mark empty strings, arrays, objects empty', function () {
				expect(angular.isEmpty('')).toBe(true);
				expect(angular.isEmpty([])).toBe(true);
				expect(angular.isEmpty({})).toBe(true);
			});

			it('should mark non-empty strings non-empty', function () {
				expect(angular.isEmpty('string')).toBe(false);
				expect(angular.isEmpty(['value'])).toBe(false);
				expect(angular.isEmpty({key : "value"})).toBe(false);
			});

			it('should work with jQuery/MooTools DOM query collections', function() {
				function Foo(elements) { Array.prototype.push.apply(this, elements); }
				Foo.prototype = { 'length': 0, 'splice': Array.prototype.splice };

				expect(angular.isEmpty(new Foo([]))).toBe(true);
			});

		});

		describe('#angular.isScope', function () {
			it('should be defined', function () {
				expect(angular.isFunction(angular.isScope)).toBe(true);
			});

			it('should mark a scope as a scope', function () {
					expect(angular.isScope($rootScope.$new())).toBe(true);
			});

			it('should mark non-scope as non-scope', function () {
				expect(angular.isScope({regular : "object"})).toBe(false);
			});
		});

		describe('#angular.once', function () {

			var count, obj;

			beforeEach(function() {
				count = 0,
					obj = {};
				obj.myOnce = angular.once(function () { count++; return count });
				spyOn(obj, 'myOnce').and.callThrough();
			});

			it('should be defined', function () {
				expect(angular.isFunction(angular.once)).toBe(true);
			});

			it('should throw a TypeError if don\'t pass function', function () {
				expect(function() { angular.once('') }).toThrow(new TypeError);
			});

			it('should return a function', function () {
				expect(angular.isFunction(obj.myOnce)).toBe(true);
			});

			it('should call the function when called', function () {
				obj.myOnce();
				expect(count).toEqual(1);
			});

			it('it should set the original function to null', function () {
				obj.myOnce();
				expect(obj.func).not.toBeDefined();
			});

			it('should only run a function once', function () {
				expect(count).toEqual(0);
				obj.myOnce();
				expect(count).toEqual(1);
				obj.myOnce();
				expect(count).toEqual(1);
			});

			it('should return same result as original function called once', function () {
				var x = obj.myOnce();
				expect(count).toEqual(1);
				expect(x).toEqual(1);
				x = obj.myOnce();
				expect(count).toEqual(1);
				expect(x).toEqual(1);
			});

			it('should maintain the `this` binding', function () {
				var func = angular.once(function() { this.count++; }),
					object = { 'count': 5, 'once': func };

				object.once();
				object.once();
				expect(object.count).toEqual(6);
			});

		});

		describe('#angular.remove', function () {
			it('should be defined', function () {
				expect(angular.isFunction(angular.remove)).toBe(true);
			});

			it('should remove objects from the original array', function () {
				var array = [1, 2, 3];

				var actual = angular.remove(array, function(num) {
					return num < 3;
				});

				expect(array).toEqual([3]);
				expect(actual).toEqual([1, 2]);
			});

			it('should support a context', function () {
				var array = [1, 2, 3];

				var actual = angular.remove(array, function(num, index) {
					return this[index] < 3;
				}, array);

				expect(actual).toEqual([1, 2]);
			});

			it('treat holes as undefined', function () {
					var array = [1, 2, 3];
					delete array[1];

					angular.remove(array, function(num) { return num == null; });
					expect(array).toEqual([1, 3]);
			});
		});

		describe('#angular.map', function () {
			it('should be defined', function () {
				expect(angular.isFunction(angular.map)).toBe(true);
			});

			it('should pass the correct arguments', function () {
				var array = [1, 2, 3];
				var args;

				angular.map(array, function() {
					args || (args = Array.prototype.slice.call(arguments));
				});

				expect(args).toEqual([1, 0, array]);
			});

			it('should support a context', function () {
				function callback(num, index) {
					return this[index] + num;
				}

				var actual = angular.map([1], callback, [2]);
				expect(actual).toEqual([3]);

				actual = angular.map({ 'a': 1 }, callback, { 'a': 2 });
				expect(actual).toEqual([3]);
			});
		});
	});
});