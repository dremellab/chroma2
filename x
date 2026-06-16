for name in scratch home standard project; do
  case "$name" in
    scratch) path="/scratch/$USER" ;;
    home) path="$HOME" ;;
    standard) path="$STANDARD" ;;
    project) path="$PROJECT" ;;
  esac

  bytes=$(
    /project/dremel_lab/cargo/bin/dust -x -o b -n 1 -c -b -P "$path" 2>/dev/null \
      | tail -n1 | awk '{print $1}' | sed 's/B$//'
  )

  printf "%s\t%s\t%s\t%s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$name" "$path" "${bytes:-NA}"
done
